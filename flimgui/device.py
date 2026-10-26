"""Device layer: a thin, defensive wrapper around pylibtdc plus a simulator.

Only the modules needed for TCSPC are touched (thresholds, channel delays,
counters, hardware histograms). The firmware updater, flash memory, oscillator
and raw I2C/DAC modules are deliberately never requested.

All device callbacks are converted to Qt signals so the GUI thread never
blocks on the USB link.
"""

from __future__ import annotations

import logging
import math
import threading
from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import QObject, QTimer, Signal

log = logging.getLogger(__name__)

SYNC_CH = 0
DEFAULT_THRESHOLD_RANGE = (0.0, 2.5)        # V, from the TDC-Software manual
DEFAULT_DELAY_RANGE = (-614e-9, 614e-9)     # s, from the TDC-Software manual
MIN_BIN_S = 10e-12                          # floor used by the simulator only


class DeviceError(RuntimeError):
    pass


@dataclass
class HistStatus:
    running: bool = False
    win_id: int = 0
    lost: bool = False
    overflow: bool = False


class DeviceBase(QObject):
    """Interface shared by the real Felix and the simulator."""

    counts_updated = Signal(float, object)           # integration time s, {channel: rate Hz}
    histogram_updated = Signal(int, object, object)  # module index, bins ndarray|None, HistStatus
    connection_lost = Signal(str)
    message = Signal(str)

    name = "device"
    is_simulator = False

    def open(self) -> None:
        raise NotImplementedError

    def close(self) -> None:
        raise NotImplementedError

    def channel_names(self) -> list[str]:
        raise NotImplementedError

    def n_channels(self) -> int:
        return len(self.channel_names())

    def n_histograms(self) -> int:
        raise NotImplementedError

    def threshold_range(self, ch: int) -> tuple[float, float]:
        return DEFAULT_THRESHOLD_RANGE

    def get_threshold(self, ch: int) -> float:
        raise NotImplementedError

    def set_threshold(self, ch: int, volts: float) -> float:
        raise NotImplementedError

    def auto_threshold(self, channels: list[int]) -> dict[int, float]:
        raise NotImplementedError

    def delay_range(self, ch: int) -> tuple[float, float]:
        return DEFAULT_DELAY_RANGE

    def get_delay(self, ch: int) -> float:
        raise NotImplementedError

    def set_delay(self, ch: int, seconds: float) -> float:
        raise NotImplementedError

    def start_counters(self, integration_s: float) -> None:
        raise NotImplementedError

    def stop_counters(self) -> None:
        raise NotImplementedError

    def hist_configure(self, idx: int, ref: int, meas: int, left_s: float, right_s: float,
                       bin_s: float | None, multi_hit: bool, accumulate: bool) -> None:
        raise NotImplementedError

    def hist_bounds(self, idx: int) -> tuple[float, float, float]:
        raise NotImplementedError

    def hist_enable(self, idx: int, refresh_s: float) -> float:
        raise NotImplementedError

    def hist_disable(self, idx: int) -> None:
        raise NotImplementedError

    def hist_running(self, idx: int) -> bool:
        raise NotImplementedError

    def hist_reset(self, idx: int) -> None:
        raise NotImplementedError

    @staticmethod
    def clamp(v: float, lo: float, hi: float) -> float:
        return max(lo, min(hi, v))


# ---------------------------------------------------------------------------
# Real hardware
# ---------------------------------------------------------------------------
class FelixDevice(DeviceBase):
    name = "Felix"

    def __init__(self, serial: str = ""):
        super().__init__()
        self.serial = serial
        self._tdcmod = None
        self._backend = None
        self._channels = []
        self._edge_requests = []
        self._thr = None
        self._counter = None
        self._hists = []
        self._tdc = None
        self._refresh = {}
        self._running = {}
        self._suppress_stop = {}
        self._lock = threading.RLock()
        self._counts_cb = None
        self._hist_cbs = {}

    # lifecycle ----------------------------------------------------------
    def open(self) -> None:
        try:
            import pylibtdc as tdc
        except ImportError as exc:
            raise DeviceError("pylibtdc is not installed in this Python") from exc
        self._tdcmod = tdc
        # NOTE: TDCBackend("", callback) accepts the arguments but aborts the process
        # (native abort, no Python frame) as soon as the device is opened, measured
        # with libTDC 1.21.8. Use the plain constructor and poll instead.
        backend = tdc.TDCBackend()
        devices = backend.getDevices()
        if not devices:
            raise DeviceError("No Felix found. Check the WinUSB driver, the USB cable and the "
                              "5 V supply, and close TDC-Software (only one program can open "
                              "the device).")
        dev = devices[0]
        if self.serial:
            for d in devices:
                if self.serial in str(list(d.values())):
                    dev = d
        err = backend.open(dev)
        if self._is_error(err):
            raise DeviceError(f"Could not open the Felix: {err}. Is TDC-Software still running?")
        backend.loadModules()
        self._backend = backend
        self._channels = list(backend.getTDCChannel())
        thr = backend.getProgrammableThreshold()
        self._thr = thr[0] if thr else None
        # The Felix has two counter modules: BELTBUS (processed events) and INPUTLOGIC
        # (events crossing the threshold at the SMA). The input rate is what the SPAD
        # really emits, so it is the one to watch for detector protection.
        ctrs = backend.getTDCCounter()
        self._counter = None
        for c in ctrs:
            if "INPUT" in str(getattr(c, "counterType", "")).upper():
                self._counter = c
        if self._counter is None and ctrs:
            self._counter = ctrs[0]
        self._hists = list(backend.getTDCHistogrammer())
        tdcs = backend.getTDC()
        self._tdc = tdcs[0] if tdcs else None
        if self._tdc is not None:
            try:
                if self._tdc.virtualDelayPresent() and not self._tdc.virtualDelayEnabled():
                    self._tdc.setVirtualDelayEnable(True)
            except Exception as exc:  # optional feature
                log.info("virtual delay enable failed: %s", exc)
        # A previous program may have left a hardware histogram running and the
        # board keeps that state. The Status returned by getAccumulator() holds
        # uninitialised values until the first acquisition of this process, so it
        # cannot tell us; stop every module unconditionally (returns at once when
        # idle, measured) and track the running state in software from here on.
        self._running = {}
        for idx, h in enumerate(self._hists):
            try:
                h.disable(True)
            except Exception as exc:
                log.warning("could not stop histogram %d at startup: %s", idx, exc)
            self._running[idx] = False
        rising = tdc.TDCChannel.edgeTrigger_type.RISING
        self._edge_requests = []
        for ch in self._channels:
            try:
                self._edge_requests.append(ch.requestEdge(rising))
            except Exception as exc:
                log.info("requestEdge failed on channel %s: %s", getattr(ch, "index", "?"), exc)
        self.message.emit(f"Connected to {dev.get('description', 'Felix')} "
                          f"S/N {dev.get('serial_number', '?')} | channels={len(self._channels)} "
                          f"histogram modules={len(self._hists)} | libTDC {backend.version()}")
        self._lost_reported = False
        self._watchdog = QTimer(self)
        self._watchdog.setInterval(2000)
        self._watchdog.timeout.connect(self._check_link)
        self._watchdog.start()

    def _check_link(self) -> None:
        """Cheap USB read; a lost link makes every call raise 'not connected'."""
        if self._backend is None or self._thr is None or self._lost_reported:
            return
        try:
            self._thr.getVoltage(0)
        except Exception as exc:
            self._lost_reported = True
            self._watchdog.stop()
            log.error("link check failed: %s", exc)
            self.connection_lost.emit(f"USB connection to the Felix was lost ({exc}).")

    @staticmethod
    def _is_error(err) -> bool:
        try:
            return bool(err.value())
        except Exception:
            try:
                return bool(err)
            except Exception:
                return False

    def close(self) -> None:
        with self._lock:
            wd = getattr(self, "_watchdog", None)
            if wd is not None:
                wd.stop()
            try:
                self.stop_counters()
            except Exception:
                pass
            for idx in range(len(self._hists)):
                try:
                    if self.hist_running(idx):
                        self._hists[idx].disable(True)
                except Exception:
                    pass
            self._edge_requests = []
            if self._backend is not None:
                try:
                    self._backend.close()
                except Exception as exc:
                    log.warning("close failed: %s", exc)
            self._backend = None

    # topology -----------------------------------------------------------
    def channel_names(self) -> list[str]:
        n = len(self._channels)
        return ["SYNC"] + [f"CH{i}" for i in range(1, n)]

    def n_histograms(self) -> int:
        return len(self._hists)

    # thresholds ---------------------------------------------------------
    def threshold_range(self, ch: int) -> tuple[float, float]:
        if self._thr is None:
            return DEFAULT_THRESHOLD_RANGE
        try:
            r = self._thr.getRange(ch)
            lo = float(getattr(r, "minVoltage", DEFAULT_THRESHOLD_RANGE[0]))
            hi = float(getattr(r, "maxVoltage", DEFAULT_THRESHOLD_RANGE[1]))
            if hi > lo:
                return (lo, hi)
        except Exception as exc:
            log.info("getRange failed: %s", exc)
        return DEFAULT_THRESHOLD_RANGE

    def get_threshold(self, ch: int) -> float:
        if self._thr is None:
            return float("nan")
        _err, v = self._thr.getVoltage(ch)
        return float(v)

    def set_threshold(self, ch: int, volts: float) -> float:
        if self._thr is None:
            raise DeviceError("no threshold module")
        lo, hi = self.threshold_range(ch)
        v = self.clamp(float(volts), lo, hi)
        err = self._thr.setVoltage(ch, v)
        if self._is_error(err):
            raise DeviceError(f"setVoltage({ch}, {v:.3f}) failed: {err}")
        return v

    def auto_threshold(self, channels: list[int]) -> dict[int, float]:
        if self._thr is None:
            raise DeviceError("no threshold module")
        out = {}
        results = self._thr.autoThreshold(list(channels))
        try:
            for item in results:
                ch, _err, v = item.get() if hasattr(item, "get") else item
                out[int(ch)] = float(v)
        except Exception as exc:
            log.info("autoThreshold result parsing: %s", exc)
        for ch in channels:
            out.setdefault(ch, self.get_threshold(ch))
        return out

    # delays -------------------------------------------------------------
    def delay_range(self, ch: int) -> tuple[float, float]:
        try:
            c = self._channels[ch]
            lo = float(c.minVirtualDelay)
            hi = float(c.maxVirtualDelay)
            if hi > lo:
                return (lo, hi)
        except Exception:
            pass
        return DEFAULT_DELAY_RANGE

    def get_delay(self, ch: int) -> float:
        try:
            return float(self._channels[ch].getVirtualDelay())
        except Exception:
            return 0.0

    def set_delay(self, ch: int, seconds: float) -> float:
        lo, hi = self.delay_range(ch)
        s = self.clamp(float(seconds), lo, hi)
        return float(self._channels[ch].setVirtualDelay(s))

    # counters -----------------------------------------------------------
    def start_counters(self, integration_s: float) -> None:
        if self._counter is None:
            raise DeviceError("no counter module")
        self._counter.setIntegrationTime(float(integration_s))
        n_ch = len(self._channels)

        def cb(t_int, counts):
            vals = list(counts)
            rates = {}
            # Measured on the Felix (3 channels): the vector has exactly one entry per
            # channel, index 0 = SYNC. Older firmware prepends an internal counter;
            # in that case the vector is one longer and we skip the first entry.
            offset = 1 if len(vals) > n_ch else 0
            for ch in range(n_ch):
                if ch + offset < len(vals) and t_int > 0:
                    rates[ch] = vals[ch + offset] / t_int
            self.counts_updated.emit(float(t_int), rates)

        self._counts_cb = cb
        self._counter.setCallBack(cb)
        self._counter.setAutoPush(True)

    def stop_counters(self) -> None:
        if self._counter is not None:
            try:
                self._counter.setAutoPush(False)
            except Exception:
                pass

    # histograms ---------------------------------------------------------
    def _h(self, idx: int):
        try:
            return self._hists[idx]
        except IndexError as exc:
            raise DeviceError(f"histogram module {idx} does not exist") from exc

    def hist_configure(self, idx, ref, meas, left_s, right_s, bin_s, multi_hit, accumulate):
        h = self._h(idx)
        if self.hist_running(idx):
            raise DeviceError("stop the histogram before changing its configuration")
        rising = self._tdcmod.TDCChannel.edgeTrigger_type.RISING
        h.setRefCh(int(ref), rising)
        h.setMeasCh(int(meas), rising)
        # Verified on the Felix (libTDC 1.21.8): the argument order is
        # (left, binWidth, right) in seconds, not (left, right, binWidth).
        h.setBounds(float(left_s), None if not bin_s else float(bin_s), float(right_s))
        h.multiHit = bool(multi_hit)
        h.accumulate = bool(accumulate)

        def cb(bins, status, _idx=idx):
            st = HistStatus(running=bool(status.running),
                            win_id=int(getattr(status, "winID", 0)),
                            lost=bool(getattr(status, "lost", False)),
                            overflow=bool(getattr(status, "overflow", False)))
            if not st.running:
                if self._suppress_stop.get(_idx):
                    # stop caused by hist_reset(): the module is re-enabled right
                    # after, so the GUI must not see a "stopped" event
                    self._suppress_stop[_idx] = False
                    return
                self._running[_idx] = False
            arr = np.array(list(bins), dtype=np.float64) if st.running else None
            self.histogram_updated.emit(_idx, arr, st)

        self._hist_cbs[idx] = cb
        h.setCallback(cb)

    def hist_bounds(self, idx: int) -> tuple[float, float, float]:
        left, bw, right = self._h(idx).getBounds()   # device order: left, binWidth, right
        return float(left), float(right), float(bw)

    def hist_enable(self, idx: int, refresh_s: float) -> float:
        h = self._h(idx)
        actual = h.enable(float(refresh_s), 0)
        self._running[idx] = True
        self._refresh[idx] = float(actual) if actual else float(refresh_s)
        return self._refresh[idx]

    def hist_disable(self, idx: int) -> None:
        self._running[idx] = False
        self._h(idx).disable(False)

    def hist_running(self, idx: int) -> bool:
        return bool(self._running.get(idx, False))

    def hist_reset(self, idx: int) -> None:
        h = self._h(idx)
        if self.hist_running(idx):
            self._suppress_stop[idx] = True
            h.disable(True)                      # blocks until the stop callback ran
            self._suppress_stop[idx] = False
            h.enable(self._refresh.get(idx, 0.5), 0)
            self._running[idx] = True


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------
def periodic_decay_s(t: np.ndarray, tau: float, period: float) -> np.ndarray:
    return np.exp(-t / tau) / (1.0 - np.exp(-period / tau))


class SimDevice(DeviceBase):
    """Fake Felix + laser + SPAD. Produces a two-exponential decay through a Gaussian IRF."""

    name = "Simulator"
    is_simulator = True
    N_MAXBIN = 8192
    LSB_S = 3.662109375e-14

    def __init__(self, n_channels: int = 3, n_hist: int = 2):
        super().__init__()
        self._n_ch = n_channels
        self._n_hist = n_hist
        self.laser_hz = 20e6
        self.sample = "fluorescence"            # or "scatter"
        self.signal_rate = {1: 3.0e5, 2: 20.0}  # Hz on each measure channel
        self.dark_rate = 20.0
        self.irf_pos_s = 8e-9
        self.irf_sigma_s = 0.12e-9
        self.components = [(0.3, 0.8e-9), (0.7, 3.6e-9)]   # (fraction, tau)
        self.pulse_amp_v = 2.5
        self._thr = {ch: 1.2 for ch in range(n_channels)}
        self._delay = {ch: 0.0 for ch in range(n_channels)}
        self._cfg = {}
        self._acc = {}
        self._win = {}
        self._timers = {}
        self._ctimer = None
        self._cint = 0.5
        self._rng = np.random.default_rng()
        self._pdf_cache = {}

    def open(self) -> None:
        self.message.emit("Simulator running: no hardware connected.")

    def close(self) -> None:
        self.stop_counters()
        for idx in list(self._timers):
            self.hist_disable(idx)

    def channel_names(self) -> list[str]:
        return ["SYNC"] + [f"CH{i}" for i in range(1, self._n_ch)]

    def n_histograms(self) -> int:
        return self._n_hist

    def get_threshold(self, ch):
        return self._thr[ch]

    def set_threshold(self, ch, volts):
        lo, hi = self.threshold_range(ch)
        self._thr[ch] = self.clamp(float(volts), lo, hi)
        return self._thr[ch]

    def auto_threshold(self, channels):
        for ch in channels:
            self._thr[ch] = self.pulse_amp_v / 2
        return {ch: self._thr[ch] for ch in channels}

    def get_delay(self, ch):
        return self._delay[ch]

    def set_delay(self, ch, seconds):
        lo, hi = self.delay_range(ch)
        self._delay[ch] = self.clamp(float(seconds), lo, hi)
        return self._delay[ch]

    def _rate(self, ch: int) -> float:
        thr = self._thr[ch]
        if thr > self.pulse_amp_v * 0.92:
            return 0.0
        noise = 2e6 * max(0.0, 0.15 - thr) / 0.15     # comparator chatter near 0 V
        if ch == SYNC_CH:
            return self.laser_hz + noise
        return self.signal_rate.get(ch, 0.0) + self.dark_rate + noise

    def start_counters(self, integration_s):
        self._cint = float(integration_s)
        if self._ctimer is None:
            self._ctimer = QTimer(self)
            self._ctimer.timeout.connect(self._tick_counts)
        self._ctimer.start(int(self._cint * 1000))

    def stop_counters(self):
        if self._ctimer is not None:
            self._ctimer.stop()

    def _tick_counts(self):
        rates = {}
        for ch in range(self._n_ch):
            mean = self._rate(ch) * self._cint
            n = self._rng.poisson(mean) if mean < 1e8 else mean
            rates[ch] = n / self._cint
        self.counts_updated.emit(self._cint, rates)

    def hist_configure(self, idx, ref, meas, left_s, right_s, bin_s, multi_hit, accumulate):
        if self.hist_running(idx):
            raise DeviceError("stop the histogram before changing its configuration")
        if right_s <= left_s:
            raise DeviceError("right bound must be larger than the left bound")
        # Mimic the Felix: bin widths are 2**k times the 36.6 fs LSB, chosen so that
        # the span fits in N_MAXBIN bins; the accumulator is always N_MAXBIN long.
        span = right_s - left_s
        k = max(math.ceil(math.log2(span / self.N_MAXBIN / self.LSB_S)), 0)
        min_bin = self.LSB_S * 2 ** k
        if not bin_s or bin_s < min_bin:
            bw = min_bin
        else:
            bw = self.LSB_S * 2 ** math.floor(math.log2(bin_s / self.LSB_S))
        n = max(int(round(span / bw)), 1)
        self._cfg[idx] = dict(ref=ref, meas=meas, left=left_s, right=left_s + n * bw, bw=bw,
                              n=n, multi_hit=multi_hit, accumulate=accumulate, refresh=0.5)
        self._acc[idx] = np.zeros(self.N_MAXBIN)
        self._win[idx] = 0

    def hist_bounds(self, idx):
        c = self._cfg[idx]
        return c["left"], c["right"], c["bw"]

    def hist_enable(self, idx, refresh_s):
        if idx not in self._cfg:
            raise DeviceError("configure the histogram first")
        self._cfg[idx]["refresh"] = max(0.001, float(refresh_s))
        if idx not in self._timers:
            t = QTimer(self)
            t.timeout.connect(lambda _i=idx: self._tick_hist(_i))
            self._timers[idx] = t
        if not self._timers[idx].isActive():
            self._acc[idx][:] = 0
            self._win[idx] = 0
        self._timers[idx].start(int(self._cfg[idx]["refresh"] * 1000))
        return self._cfg[idx]["refresh"]

    def hist_disable(self, idx):
        if self.hist_running(idx):
            self._timers[idx].stop()
            self.histogram_updated.emit(idx, None, HistStatus(running=False))

    def hist_running(self, idx):
        return idx in self._timers and self._timers[idx].isActive()

    def hist_reset(self, idx):
        if self.hist_running(idx):
            self._acc[idx][:] = 0
            self._win[idx] = 0

    def _period_pdf(self):
        """Probability density over one laser period on a fine grid (cached)."""
        T = 1.0 / self.laser_hz
        key = (self.sample, T, self.irf_sigma_s, tuple(self.components))
        if key in self._pdf_cache:
            return self._pdf_cache[key]
        m = 4096
        tt = (np.arange(m) + 0.5) * T / m
        irf = np.exp(-0.5 * ((tt - 4 * self.irf_sigma_s) / self.irf_sigma_s) ** 2)
        irf /= irf.sum()
        if self.sample == "scatter":
            pdf = irf.copy()
        else:
            dec = np.zeros(m)
            for frac, tau in self.components:
                dec += frac / tau * periodic_decay_s(tt, tau, T)
            pdf = np.fft.irfft(np.fft.rfft(irf) * np.fft.rfft(dec), n=m)
            pdf = np.clip(pdf, 0, None)
        pdf /= pdf.sum() * (T / m)       # density per second over one period
        self._pdf_cache[key] = (tt, pdf, T)
        return self._pdf_cache[key]

    def _tick_hist(self, idx):
        c = self._cfg[idx]
        meas, ref = c["meas"], c["ref"]
        n, bw, left = c["n"], c["bw"], c["left"]
        tt, pdf, T = self._period_pdf()
        centers = left + (np.arange(n) + 0.5) * bw
        offset = self.irf_pos_s + self._delay.get(meas, 0.0) - self._delay.get(ref, 0.0)
        phase = np.mod(centers - offset, T)
        dens = np.interp(phase, tt, pdf, period=T)
        sig_rate = self.signal_rate.get(meas, 0.0) if self._rate(meas) > 0 else 0.0
        expected = sig_rate * c["refresh"] * dens * bw
        bg_per_bin = (self.dark_rate + 0.002 * sig_rate) * c["refresh"] * bw / T
        expected = expected + bg_per_bin
        window = np.zeros(self.N_MAXBIN)
        window[:n] = self._rng.poisson(np.clip(expected, 0, None))
        if c["accumulate"]:
            self._acc[idx] += window
            out = self._acc[idx]
        else:
            out = window
        self._win[idx] += 1
        st = HistStatus(running=True, win_id=self._win[idx], lost=False,
                        overflow=bool(out.max() >= 2 ** 32 - 1))
        self.histogram_updated.emit(idx, out.copy(), st)
