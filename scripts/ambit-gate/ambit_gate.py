"""AMBIT + Felix: fluorescence yield and software-gated lifetimes around a saturating flash.

What it does
------------
1. Finds the AMBIT on a serial port (text console, ``hello``).
2. Finds the Felix TDC with pylibtdc, sets the input thresholds, prints count rates.
3. Starts the Felix hardware histogram in period mode (contiguous windows of
   ``window_s``, e.g. 100 ms), waits ``pre_s``, sends the ``arrun`` protocol to the
   AMBIT and collects both data sets.
4. Aligns the two clocks. The SPAD saturates while the actinic LED is on, so the
   window totals show the flash on/off edges directly. The AMBIT trace is placed
   relative to those edges; the nominal command latency is the fallback.
5. Classifies every Felix window (pre, flash, post, edge, outside), computes a
   mean lifetime per window, and writes everything to an output folder with a plot.

"Software gating" means there is no cable between the AMBIT and the Felix: the
windows that overlap the flash are discarded afterwards. See README.md.

Conventions verified on the real Felix (see flimgui/device.py):
  * setBounds argument order is (left, binWidth, right) in seconds;
  * the counter vector has one entry per channel, index 0 = SYNC;
  * the callback always delivers 8192 bins, only the first n are real.

Run ``python ambit_gate.py --check`` to detect both instruments without measuring.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import serial
import serial.tools.list_ports

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))  # for flimgui.fitting (optional)

try:
    from flimgui.fitting import fit_decay  # tail fit, no IRF needed
except Exception:  # scipy missing or repo moved
    fit_decay = None

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BAUD = 115200
LASER_REP_RATE_HZ = 80e6  # laser mode selector; 12.5 ns period
SYNC_CH, SPAD_CH = 0, 1  # Felix SYNC = channel 0, SPAD on CH1
THRESHOLD_V = 1.2  # both signals are 2.5 V LVTTL
WINDOW_S = 0.1  # Felix window length (period mode)
PRE_S = 1.0  # Felix runs this long before the arrun is sent
POST_S = 1.0  # ... and this long after the protocol ends
NOMINAL_LATENCY_S = 0.04  # command parse + config + warm-up before sample 0 (fallback)
ACTINIC_ON = 3  # firmware: actinic byte > 3 switches the LED on

# The experiment asked for: 1 s dark, 1 s actinic at 250 mA, 3 s relaxation, 10 Hz.
# type 2 = no 730 nm reference pulses (keeps the SPAD band clean).
DEFAULT_PROTOCOL = [
    dict(type=2, far_red=0, pulses=10, freq=10, actinic=0, subsampling=0),
    dict(type=2, far_red=0, pulses=10, freq=10, actinic=250, subsampling=0),
    dict(type=2, far_red=0, pulses=30, freq=10, actinic=0, subsampling=0),
]

# USB-serial bridges an AMBIT is usually seen behind; probed first.
KNOWN_VIDS = {0x1A86: "WCH CH340/CH9102", 0x10C4: "Silicon Labs CP210x", 0x0403: "FTDI", 0x303A: "Espressif native USB"}
HELLO_RE = re.compile(r"NEW\s+(?P<name>.*?)\s+Ready(?:\s+FW:(?P<fw>\S+))?")
DATA_RE = re.compile(r"Data:(?P<tag>[^,]+),Length:(?P<n>\d+)\s*(?P<vals>.*)")
MAX_SEGMENTS = 16  # do_command.h: arr[128] = 16 lines x 8 bytes
MAX_POINTS = 1999  # data_utils.h: MAX_DATACLASS_SIZE - 1


# ---------------------------------------------------------------------------
# AMBIT (text console over UART0)
# ---------------------------------------------------------------------------
@dataclass
class AmbitLink:
    ser: serial.Serial
    name: str
    fw: str

    @property
    def port(self) -> str:
        return self.ser.port

    def close(self) -> None:
        try:
            self.ser.close()
        except Exception:
            pass


def _open_port(device: str, timeout: float = 0.2) -> serial.Serial:
    """Open without toggling DTR/RTS so the ESP32-C3 is not reset into boot mode."""
    ser = serial.Serial()
    ser.port = device
    ser.baudrate = BAUD
    ser.timeout = timeout
    ser.dtr = False
    ser.rts = False
    ser.open()
    return ser


def _say_hello(ser: serial.Serial, wait: float = 1.5):
    ser.reset_input_buffer()
    ser.write(b"hello\n")
    deadline = time.monotonic() + wait
    while time.monotonic() < deadline:
        line = ser.readline().decode("ascii", errors="replace").strip()
        m = HELLO_RE.search(line)
        if m:
            return m
    return None


def _probe(device: str, log=print):
    try:
        ser = _open_port(device)
    except serial.SerialException as exc:
        log(f"  {device}: cannot open ({exc})")
        return None
    time.sleep(0.3)  # drop a boot banner if opening reset the board
    m = _say_hello(ser)
    if m is None:
        time.sleep(2.0)  # give a reset board time to boot, ask again
        m = _say_hello(ser)
    if m is None:
        ser.close()
        log(f"  {device}: no AMBIT reply")
        return None
    return AmbitLink(ser=ser, name=m["name"], fw=m["fw"] or "")


def find_ambit(port: str | None = None, log=print) -> AmbitLink:
    """Probe serial ports with ``hello`` until an AMBIT answers."""
    if port:
        candidates = [port]
    else:
        ports = sorted(serial.tools.list_ports.comports(), key=lambda p: p.vid not in KNOWN_VIDS)
        for p in ports:
            log(f"  found {p.device}: {p.description}")
        candidates = [p.device for p in ports]
    for dev in candidates:
        link = _probe(dev, log)
        if link:
            log(f"AMBIT '{link.name}' on {dev}, firmware {link.fw or '(unreported)'}")
            return link
    raise RuntimeError("No AMBIT found. Check the cable, or pass the port explicitly.")


def build_arrun(protocol: list[dict], persist: int = 0, verb: str = "arrun2") -> str:
    """Validate the segments and return the console command string.

    ``arrun2`` selects the COMPUTER sink: nothing is printed during the run and the
    arrays are dumped at the end, so the serial line is silent while measuring.
    """
    if not 1 <= len(protocol) <= MAX_SEGMENTS:
        raise ValueError(f"need 1..{MAX_SEGMENTS} segments, got {len(protocol)}")
    if sum(s["pulses"] for s in protocol) > MAX_POINTS:
        raise ValueError(f"total pulses must be <= {MAX_POINTS}")
    values = [len(protocol), int(bool(persist))]
    for i, s in enumerate(protocol):
        if s["type"] not in (1, 2):
            raise ValueError(f"segment {i}: type must be 1 or 2")
        if not 1 <= s["pulses"] <= 0xFFFF or not 1 <= s["freq"] <= 2000:
            raise ValueError(f"segment {i}: pulses 1..65535, freq 1..2000 Hz")
        if not 0 <= s["actinic"] <= 255 or s["subsampling"] not in (0, 1, 2):
            raise ValueError(f"segment {i}: actinic 0..255, subsampling 0/1/2")
        values += [
            s["type"],
            s["far_red"],
            s["pulses"] >> 8,
            s["pulses"] & 0xFF,
            s["freq"] >> 8,
            s["freq"] & 0xFF,
            s["actinic"],
            s["subsampling"],
        ]
    # Every argument is read with a ',' terminator, so the line ends with a comma.
    return f"{verb}," + ",".join(map(str, values)) + ",\n"


def protocol_duration(protocol: list[dict]) -> float:
    return sum(s["pulses"] / s["freq"] for s in protocol)


def time_axis(protocol: list[dict]):
    """Nominal seconds since sample 0 for every fluorescence point, plus segment edges.

    The firmware paces each edge from the previous one with a period of 1/freq, so
    the nominal axis is accurate to microseconds; only the run start is uncertain.
    """
    t, edges, t0 = [], [], 0.0
    for s in protocol:
        n, f = s["pulses"], s["freq"]
        t.append(t0 + np.arange(n) / f)
        edges.append((t0, t0 + n / f, s))
        t0 += n / f
    return np.concatenate(t), edges


def run_arrun(link: AmbitLink, cmd: str, protocol: list[dict], margin_s: float = 10.0, log=print):
    """Send ``cmd``, read until ``Data sent``. Returns (data dict, t_send, t_done)."""
    duration = protocol_duration(protocol)
    ser = link.ser
    ser.reset_input_buffer()
    t_send = time.monotonic()
    ser.write(cmd.encode("ascii"))
    ser.flush()
    deadline = t_send + duration + margin_s
    data, other = {}, []
    while time.monotonic() < deadline:
        line = ser.readline().decode("ascii", errors="replace").strip()
        if not line:
            continue
        if line == "Data sent":
            break
        if line.startswith("ERROR"):
            raise RuntimeError(f"the AMBIT refused the run: {line}")
        if line.startswith("ESP-ROM") or line.startswith("rst:"):
            raise RuntimeError("the AMBIT rebooted during the run")
        m = DATA_RE.match(line)
        if m:
            vals = [int(v) for v in m["vals"].split(",") if v.strip()]
            if len(vals) != int(m["n"]):
                log(f"warning: {m['tag']} announced {m['n']} values, got {len(vals)}")
            data[m["tag"]] = np.array(vals, dtype=np.int64)
        else:
            other.append(line)
    else:
        raise TimeoutError(f"no 'Data sent' after {duration + margin_s:.0f} s; got {list(data)}")
    t_done = time.monotonic()
    if other:
        log("other console output: " + " | ".join(other))
    return data, t_send, t_done


# ---------------------------------------------------------------------------
# Felix (pylibtdc, no Qt)
# ---------------------------------------------------------------------------
@dataclass
class Window:
    k: int  # window index (0 = first after start of acquisition)
    win_id: int
    t_cb: float  # monotonic time the window arrived on the PC
    bins: np.ndarray  # real bins only (uint64)
    lost: bool = False
    overflow: bool = False


FELIX_VID, FELIX_PID = 0x0403, 0x6C68  # USB ids of the Felix (TEDIEL); interface 0 is the WinUSB one
CM_PROB_FAILED_INSTALL = 28  # Windows problem code: device present, no driver installed


def felix_usb_state(interface: int = 0) -> dict:
    """What Windows itself knows about the Felix USB interface, independent of pylibtdc.

    Returns ``{"present": bool | None, "problem": int | None, "instance": str}``: ``present`` is
    None on other platforms or when the device tree cannot be read, ``problem`` is the Device
    Manager problem code (0 = fine, 28 = no driver installed). Lets "No Felix found" say why.
    """
    none = {"present": None, "problem": None, "instance": ""}
    if sys.platform != "win32":
        return none
    import ctypes
    import ctypes.wintypes as wt

    class DevInfoData(ctypes.Structure):
        _fields_ = [
            ("cbSize", wt.DWORD),
            ("ClassGuid", ctypes.c_byte * 16),
            ("DevInst", wt.DWORD),
            ("Reserved", ctypes.c_void_p),
        ]

    setupapi, cfgmgr = ctypes.WinDLL("setupapi"), ctypes.WinDLL("cfgmgr32")
    setupapi.SetupDiGetClassDevsW.restype = ctypes.c_void_p
    setupapi.SetupDiGetClassDevsW.argtypes = [ctypes.c_void_p, wt.LPCWSTR, wt.HWND, wt.DWORD]
    setupapi.SetupDiEnumDeviceInfo.argtypes = [ctypes.c_void_p, wt.DWORD, ctypes.POINTER(DevInfoData)]
    setupapi.SetupDiGetDeviceInstanceIdW.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(DevInfoData),
        wt.LPWSTR,
        wt.DWORD,
        ctypes.POINTER(wt.DWORD),
    ]
    setupapi.SetupDiDestroyDeviceInfoList.argtypes = [ctypes.c_void_p]
    cfgmgr.CM_Get_DevNode_Status.argtypes = [ctypes.POINTER(wt.ULONG), ctypes.POINTER(wt.ULONG), wt.DWORD, wt.ULONG]
    DIGCF_PRESENT, DIGCF_ALLCLASSES = 0x2, 0x4
    hdev = setupapi.SetupDiGetClassDevsW(None, "USB", None, DIGCF_PRESENT | DIGCF_ALLCLASSES)
    if not hdev or hdev == ctypes.c_void_p(-1).value:
        return none
    needle = f"VID_{FELIX_VID:04X}&PID_{FELIX_PID:04X}&MI_{interface:02d}"
    try:
        i = 0
        while True:
            dd = DevInfoData(cbSize=ctypes.sizeof(DevInfoData))
            if not setupapi.SetupDiEnumDeviceInfo(hdev, i, ctypes.byref(dd)):
                return {"present": False, "problem": None, "instance": ""}
            i += 1
            buf, n = ctypes.create_unicode_buffer(512), wt.DWORD()
            if not setupapi.SetupDiGetDeviceInstanceIdW(hdev, ctypes.byref(dd), buf, 512, ctypes.byref(n)):
                continue
            if needle not in buf.value.upper():
                continue
            status, problem = wt.ULONG(), wt.ULONG()
            cfgmgr.CM_Get_DevNode_Status(ctypes.byref(status), ctypes.byref(problem), dd.DevInst, 0)
            return {"present": True, "problem": int(problem.value), "instance": buf.value}
    finally:
        setupapi.SetupDiDestroyDeviceInfoList(hdev)


def felix_usb_hint() -> str:
    """One sentence on why pylibtdc sees no Felix, from the Windows device tree."""
    st = felix_usb_state()
    if st["present"] is None:
        return "Check the USB cable and the 5 V supply, and close TDC-Software (one program at a time)."
    if not st["present"]:
        return "No Felix on USB: check the USB cable and the 5 V supply."
    if st["problem"] == CM_PROB_FAILED_INSTALL:
        return (
            "The Felix is plugged in but Windows has no driver for it: install TEDIEL's Felix WinUSB "
            "driver (the phi2 GUI offers to do it; also at https://tediel.com/resources)."
        )
    if st["problem"]:
        return f"The Felix is plugged in but Windows reports problem code {st['problem']} for it (Device Manager)."
    return (
        "The Felix is plugged in with its driver but could not be opened: close TDC-Software "
        "(only one program can use the Felix), or unplug and replug it."
    )


REJECTION_WARN = 0.5  # TDC drop fraction above which the Felix is declared wedged


def rejection_warning(rejection: float, spad_cps: float = float("inf")) -> str | None:
    """Warning text when ``FelixLink.tdc_rejection`` says the TDC drops the SPAD events, else None
    (ignored below 10 kcps, where the fraction is noise)."""
    if rejection != rejection or rejection <= REJECTION_WARN or spad_cps < 1e4:
        return None
    return (
        f"the Felix TDC drops {rejection:.0%} of the SPAD events: the histograms would be empty. "
        "Unplug and replug the Felix (or power-cycle it) and reconnect"
    )


class FelixLink:
    def __init__(self, serial_hint: str = "", log=print):
        self.log = log
        self.serial_hint = serial_hint
        self.tdc = None
        self.backend = None
        self.channels = []
        self._edge_requests = []
        self.thr = None
        self.counter = None
        self.hist = None
        self.n_bins = 0
        self.left_s = 0.0
        self.right_s = 0.0
        self.bin_s = 0.0
        self.windows: list[Window] = []
        self._lock = threading.Lock()
        self._done = threading.Event()
        self._n_expected = 0
        self.t_soa = float("nan")
        self.window_s = float("nan")
        self.description = ""

    @staticmethod
    def _is_error(err) -> bool:
        try:
            return bool(err.value())
        except Exception:
            try:
                return bool(err)
            except Exception:
                return False

    def open(self) -> FelixLink:
        try:
            import pylibtdc as tdc
        except ImportError as exc:
            raise RuntimeError(
                "pylibtdc is not installed in this Python: run 'uv sync' in the repository "
                "folder. TEDIEL ships it for Windows and Linux (x64) only, so a Mac cannot "
                "drive the Felix"
            ) from exc
        self.tdc = tdc
        backend = tdc.TDCBackend()  # never pass a callback: it aborts the process
        devices = backend.getDevices()
        if not devices:
            raise RuntimeError("No Felix found. " + felix_usb_hint())
        dev = devices[0]
        if self.serial_hint:
            for d in devices:
                if self.serial_hint in str(list(d.values())):
                    dev = d
        err = backend.open(dev)
        if self._is_error(err):
            raise RuntimeError(f"Could not open the Felix: {err}. Is TDC-Software running?")
        backend.loadModules()
        self.backend = backend
        self.channels = list(backend.getTDCChannel())
        thr = backend.getProgrammableThreshold()
        self.thr = thr[0] if thr else None
        self.counter = None
        for c in backend.getTDCCounter():  # INPUTLOGIC = rate at the SMA
            if "INPUT" in str(getattr(c, "counterType", "")).upper():
                self.counter = c
        if self.counter is None:
            ctrs = backend.getTDCCounter()
            self.counter = ctrs[0] if ctrs else None
        hists = list(backend.getTDCHistogrammer())
        if not hists:
            raise RuntimeError("the Felix reports no histogram module")
        for h in hists:  # a previous program may have left one running
            try:
                h.disable(True)
            except Exception:
                pass
        self.hist = hists[0]
        tdcs = backend.getTDC()
        if tdcs:
            try:
                if tdcs[0].virtualDelayPresent() and not tdcs[0].virtualDelayEnabled():
                    tdcs[0].setVirtualDelayEnable(True)
            except Exception:
                pass
        rising = tdc.TDCChannel.edgeTrigger_type.RISING
        self._edge_requests = [ch.requestEdge(rising) for ch in self.channels]
        self.description = (
            f"{dev.get('description', 'Felix')} S/N {dev.get('serial_number', '?')}"
            f" | {len(self.channels)} channels, {len(hists)} histogram modules"
            f" | libTDC {backend.version()}"
        )
        self.log("Felix: " + self.description)
        return self

    def close(self) -> None:
        try:
            if self.hist is not None:
                self.hist.disable(True)
        except Exception:
            pass
        self._edge_requests = []
        if self.backend is not None:
            try:
                self.backend.close()
            except Exception:
                pass
        self.backend = None

    # configuration --------------------------------------------------------
    def set_thresholds(self, volts: float, channels=(SYNC_CH, SPAD_CH)) -> None:
        if self.thr is None:
            raise RuntimeError("no threshold module")
        for ch in channels:
            err = self.thr.setVoltage(ch, float(volts))
            if self._is_error(err):
                raise RuntimeError(f"setVoltage({ch}, {volts}) failed: {err}")

    def count_rates(self, integration_s: float = 0.5) -> dict[int, float]:
        """Input rate (Hz) per channel, index 0 = SYNC."""
        if self.counter is None:
            raise RuntimeError("no counter module")
        self.counter.setIntegrationTime(float(integration_s))
        t_int, counts = self.counter.getCounts()
        vals = list(counts)
        n_ch = len(self.channels)
        offset = 1 if len(vals) > n_ch else 0  # older firmware prepends an internal counter
        return {ch: vals[ch + offset] / t_int for ch in range(n_ch) if ch + offset < len(vals) and t_int > 0}

    def tdc_rejection(self, integration_s: float = 0.5, channel: int = SPAD_CH) -> float:
        """Fraction of the events arriving on ``channel`` that the TDC core drops: growth of the
        channel's dead-time counter against the input counter over ``integration_s``.

        About 0 when healthy. 1.0 is the wedged state seen on 2026-11-05: input rates looked fine
        while every histogram window came back empty ("no laser peak"); a power cycle of the Felix
        cleared it, no register had changed. NaN when the registers cannot be read.
        """
        try:
            ch = self.channels[channel]
            prop = ch.channel_property
            d0 = ch.get(prop.DEADCOUNT_R)
            self.counter.setIntegrationTime(float(integration_s))
            _t_int, counts = self.counter.getCounts()
            d1 = ch.get(prop.DEADCOUNT_R)
            vals = list(counts)
            offset = 1 if len(vals) > len(self.channels) else 0
            n_in = vals[channel + offset] if channel + offset < len(vals) else 0
            return ((d1 - d0) % (1 << 32)) / n_in if n_in > 1000 else 0.0
        except Exception:
            return float("nan")

    def configure_histogram(
        self, ref=SYNC_CH, meas=SPAD_CH, left_s=0.0, right_s=1.0 / LASER_REP_RATE_HZ, bin_s=None
    ) -> None:
        rising = self.tdc.TDCChannel.edgeTrigger_type.RISING
        h = self.hist
        h.setRefCh(int(ref), rising)
        h.setMeasCh(int(meas), rising)
        h.setBounds(float(left_s), None if not bin_s else float(bin_s), float(right_s))
        h.multiHit = False
        h.accumulate = False  # every window is delivered on its own
        left, bw, right = h.getBounds()  # device order: left, binWidth, right
        self.left_s, self.bin_s, self.right_s = float(left), float(bw), float(right)
        self.n_bins = int(round((self.right_s - self.left_s) / self.bin_s))
        self.log(
            f"Felix histogram: {self.left_s * 1e9:.2f}..{self.right_s * 1e9:.2f} ns, "
            f"bin {self.bin_s * 1e12:.2f} ps, {self.n_bins} bins"
        )

    # acquisition ----------------------------------------------------------
    def start_windows(self, window_s: float, n_windows: int) -> float:
        """Start period-mode acquisition: ``n_windows`` contiguous windows of ``window_s``.

        Returns the actual window length. ``self.t_soa`` is the PC monotonic time of
        the start of acquisition (enable() waits for it).
        """
        self.windows = []
        self._done.clear()
        self._n_expected = int(n_windows)
        n_bins = self.n_bins

        def cb(bins, status):
            running = bool(getattr(status, "running", True))
            if not running:
                self._done.set()
                return
            arr = np.array(list(bins)[:n_bins], dtype=np.uint64)
            with self._lock:
                w = Window(
                    k=len(self.windows),
                    win_id=int(getattr(status, "winID", 0)),
                    t_cb=time.monotonic(),
                    bins=arr,
                    lost=bool(getattr(status, "lost", False)),
                    overflow=bool(getattr(status, "overflow", False)),
                )
                self.windows.append(w)
                if len(self.windows) >= self._n_expected:
                    self._done.set()

        self.hist.setCallback(cb)
        try:
            self.hist.enableWaitForSOA(2.0)  # enable() returns at the real start
        except Exception:
            pass
        actual = self.hist.enable(float(window_s), int(n_windows))
        self.t_soa = time.monotonic()
        self.window_s = float(actual) if actual else float(window_s)
        return self.window_s

    def wait_done(self, timeout_s: float) -> bool:
        ok = self._done.wait(timeout_s)
        try:
            self.hist.disable(True)
        except Exception:
            pass
        try:
            self.hist.setCallback()  # clear the user callback
        except Exception:
            pass
        return ok


# ---------------------------------------------------------------------------
# Measurement, alignment, analysis
# ---------------------------------------------------------------------------
@dataclass
class RunResult:
    protocol: list[dict]
    cmd: str
    ambit: dict  # tag -> array
    windows: list[Window]
    window_s: float
    t_soa: float
    t_send: float
    t_done: float
    bin_s: float
    left_s: float
    laser_hz: float
    rates_before: dict
    ambit_name: str = ""
    ambit_fw: str = ""
    felix_desc: str = ""
    run_start: float = float("nan")  # Felix-clock (monotonic) time of AMBIT sample 0
    align_mode: str = "nominal"
    flash_measured: tuple = ()  # (onset, end) in monotonic seconds, if detected
    summary: list[dict] = field(default_factory=list)

    # derived ---------------------------------------------------------------
    @property
    def win_start(self) -> np.ndarray:
        return self.t_soa + np.array([w.k for w in self.windows]) * self.window_s

    @property
    def totals(self) -> np.ndarray:
        return np.array([w.bins.sum() for w in self.windows], dtype=float)

    def ambit_time(self) -> np.ndarray:
        """Monotonic time of every AMBIT fluorescence sample."""
        t, _ = time_axis(self.protocol)
        return self.run_start + t

    def bins_matrix(self) -> np.ndarray:
        return np.vstack([w.bins for w in self.windows]) if self.windows else np.zeros((0, 0))


def run_gated(
    ambit: AmbitLink,
    felix: FelixLink,
    protocol=None,
    window_s=WINDOW_S,
    pre_s=PRE_S,
    post_s=POST_S,
    laser_hz=LASER_REP_RATE_HZ,
    threshold_v=THRESHOLD_V,
    out_dir: Path | str | None = None,
    log=print,
) -> RunResult:
    """Run the protocol on the AMBIT while the Felix records contiguous windows."""
    protocol = list(protocol or DEFAULT_PROTOCOL)
    cmd = build_arrun(protocol)
    duration = protocol_duration(protocol)
    n_windows = int(np.ceil((pre_s + duration + post_s) / window_s))

    felix.set_thresholds(threshold_v)
    rates = felix.count_rates(0.5)
    log("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in rates.items()))
    if rates.get(SYNC_CH, 0) < 0.9 * laser_hz:
        log("WARNING: SYNC rate below the laser repetition rate. Check SYNC cable, laser mode, threshold.")
    warn = rejection_warning(felix.tdc_rejection(0.5), rates.get(SPAD_CH, 0))
    if warn:
        log("WARNING: " + warn)
    felix.configure_histogram(SYNC_CH, SPAD_CH, 0.0, 1.0 / laser_hz)

    log(f"Felix: {n_windows} windows of {window_s * 1e3:.0f} ms ...")
    actual = felix.start_windows(window_s, n_windows)
    if abs(actual - window_s) > 0.05 * window_s:
        log(f"note: Felix rounded the window to {actual * 1e3:.2f} ms")
    time.sleep(pre_s)
    log("AMBIT: " + cmd.strip())
    data, t_send, t_done = run_arrun(ambit, cmd, protocol, log=log)
    ok = felix.wait_done(timeout_s=post_s + 5.0)
    if not ok:
        log(f"WARNING: Felix delivered {len(felix.windows)} of {n_windows} windows")

    res = RunResult(
        protocol=protocol,
        cmd=cmd,
        ambit=data,
        windows=list(felix.windows),
        window_s=felix.window_s,
        t_soa=felix.t_soa,
        t_send=t_send,
        t_done=t_done,
        bin_s=felix.bin_s,
        left_s=felix.left_s,
        laser_hz=laser_hz,
        rates_before=rates,
        ambit_name=ambit.name,
        ambit_fw=ambit.fw,
        felix_desc=felix.description,
    )
    align(res, log=log)
    analyse(res, log=log)
    if out_dir is not None:
        save_run(res, Path(out_dir), log=log)
    return res


def align(res: RunResult, nominal_latency_s: float = NOMINAL_LATENCY_S, log=print) -> None:
    """Place the AMBIT time axis on the Felix clock.

    Preferred: the flash saturates the SPAD, so the first and last windows with a
    high total count mark the flash, and the partially filled neighbours refine the
    edge to a few ms. Fallback: command send time plus the nominal latency.
    """
    res.run_start = res.t_send + nominal_latency_s
    res.align_mode = "nominal"
    _, edges = time_axis(res.protocol)
    flash = [(t0, t1) for t0, t1, s in edges if s["actinic"] > ACTINIC_ON]
    if not flash or len(res.windows) < 5:
        return
    tot = res.totals
    starts = res.win_start
    W = res.window_s
    if any(w.lost for w in res.windows):
        log("WARNING: the Felix reported lost windows; window timing is unreliable")
    # Fully lit windows define the saturated level, windows below half of it the
    # dark level. The first and last windows above half level may be partially
    # lit; their fill fraction places the edge inside the window.
    sat = float(np.median(tot[tot > 0.9 * tot.max()]))
    lit = tot > 0.5 * sat
    if lit.sum() == 0 or (~lit).sum() == 0:
        return
    dark = float(np.median(tot[~lit]))
    if dark <= 0 or sat / dark < 3.0:
        log(
            f"flash not visible in the Felix totals (contrast {sat / max(dark, 1):.1f}); "
            f"using the nominal latency of {nominal_latency_s * 1e3:.0f} ms"
        )
        return
    i0 = int(np.argmax(lit))  # first window with light
    i1 = i0
    while i1 + 1 < len(tot) and lit[i1 + 1]:
        i1 += 1
    frac = lambda c: float(np.clip((c - dark) / (sat - dark), 0.0, 1.0))
    onset = starts[i0] + (1.0 - frac(tot[i0])) * W
    end = starts[i1] + frac(tot[i1]) * W
    f0, f1 = flash[0]
    res.run_start = onset - f0
    res.align_mode = "flash-edges"
    res.flash_measured = (onset, end)
    latency = res.run_start - res.t_send
    log(
        f"alignment from the flash edges: run start = send + {latency * 1e3:.0f} ms, "
        f"flash lasted {end - onset:.3f} s (protocol {f1 - f0:.3f} s)"
    )
    if abs((end - onset) - (f1 - f0)) > 0.15:
        log("WARNING: measured flash length differs from the protocol; check the alignment")
    if not (0.0 <= latency <= 0.5):
        log("WARNING: implausible latency; falling back to the nominal value")
        res.run_start = res.t_send + nominal_latency_s
        res.align_mode = "nominal"


def classify_window(t0: float, t1: float, run_start: float, edges) -> str:
    """pre / flash / post / edge / outside, for a window [t0, t1] on the Felix clock."""
    seg_hits = []
    for i, (a, b, s) in enumerate(edges):
        a, b = run_start + a, run_start + b
        if t1 <= a or t0 >= b:
            continue
        seg_hits.append((i, s, t0 >= a and t1 <= b))
    if not seg_hits:
        return "outside"
    if len(seg_hits) > 1 or not seg_hits[0][2]:
        return "edge"
    i, s, _ = seg_hits[0]
    if s["actinic"] > ACTINIC_ON:
        return "flash"
    first_flash = next((j for j, (_, _, ss) in enumerate(edges) if ss["actinic"] > ACTINIC_ON), None)
    if first_flash is None or i < first_flash:
        return "pre"
    return "post"


def mean_lifetime_ns(bins: np.ndarray, bin_s: float) -> tuple[float, float]:
    """(mean arrival time after the peak, background per bin) with a flat background removed.

    A model-free proxy for the lifetime; use ``fit_decay`` for a real fit.
    """
    c = np.asarray(bins, float)
    if c.sum() < 50:
        return float("nan"), float("nan")
    bg = float(np.percentile(c, 5))
    peak = int(np.argmax(c))
    tail = c[peak:] - bg
    t = np.arange(tail.size) * bin_s * 1e9
    s = tail.sum()
    return (float((tail * t).sum() / s) if s > 0 else float("nan")), bg


def analyse(res: RunResult, n_exp: int = 1, log=print) -> list[dict]:
    _, edges = time_axis(res.protocol)
    rows = []
    period_ns = 1e9 / res.laser_hz
    for w, t0 in zip(res.windows, res.win_start):
        t1 = t0 + res.window_s
        phase = classify_window(t0, t1, res.run_start, edges)
        tau_mean, bg = mean_lifetime_ns(w.bins, res.bin_s)
        tau_fit = float("nan")
        if fit_decay is not None and phase in ("pre", "post") and w.bins.sum() >= 200:
            try:
                fr = fit_decay(w.bins.astype(float), res.bin_s * 1e9, n_exp=n_exp, period_ns=period_ns)
                if fr.success:
                    tau_fit = fr.tau_int_ns
            except Exception as exc:
                log(f"fit failed on window {w.k}: {exc}")
        rows.append(
            dict(
                k=w.k,
                win_id=w.win_id,
                t_start_s=t0 - res.run_start,
                t_end_s=t1 - res.run_start,
                phase=phase,
                total=int(w.bins.sum()),
                bg_per_bin=bg,
                tau_mean_ns=tau_mean,
                tau_fit_ns=tau_fit,
                lost=w.lost,
                overflow=w.overflow,
            )
        )
    res.summary = rows
    n = {p: sum(r["phase"] == p for r in rows) for p in ("pre", "flash", "post", "edge", "outside")}
    log("windows: " + ", ".join(f"{k}={v}" for k, v in n.items()))
    return rows


def summed_decay(res: RunResult, phase: str) -> np.ndarray:
    sel = [w.bins for w, r in zip(res.windows, res.summary) if r["phase"] == phase]
    if not sel:
        return np.zeros(res.windows[0].bins.size if res.windows else 0, dtype=np.uint64)
    return np.sum(np.vstack(sel), axis=0)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def save_run(res: RunResult, out_dir: Path, log=print) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t_amb = res.ambit_time() - res.run_start
    s = res.ambit.get("s_630", np.zeros(0))
    r = res.ambit.get("r_630", np.zeros(0))
    n = min(len(t_amb), len(s), len(r))
    _, edges = time_axis(res.protocol)
    with open(out_dir / "ambit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["t_s", "s_630", "r_630", "ratio", "actinic"])
        for i in range(n):
            act = next((seg["actinic"] for a, b, seg in edges if a <= t_amb[i] < b), 0)
            ratio = s[i] / r[i] if r[i] > 0 else float("nan")
            wr.writerow([f"{t_amb[i]:.4f}", s[i], r[i], f"{ratio:.5f}", act])
    with open(out_dir / "felix_windows.csv", "w", newline="") as fh:
        keys = list(res.summary[0].keys()) if res.summary else []
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(res.summary)
    np.savez_compressed(
        out_dir / "felix_windows.npz",
        bins=res.bins_matrix(),
        t_start_s=np.array([r_["t_start_s"] for r_ in res.summary]),
        phase=np.array([r_["phase"] for r_ in res.summary]),
        bin_s=res.bin_s,
        left_s=res.left_s,
        window_s=res.window_s,
    )
    for phase in ("pre", "post"):
        dec = summed_decay(res, phase)
        with open(out_dir / f"decay_{phase}.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["t_ns", "counts"])
            for i, c in enumerate(dec):
                wr.writerow([f"{(res.left_s + (i + 0.5) * res.bin_s) * 1e9:.4f}", int(c)])
    meta = dict(
        protocol=res.protocol,
        cmd=res.cmd.strip(),
        window_s=res.window_s,
        laser_hz=res.laser_hz,
        bin_s=res.bin_s,
        align_mode=res.align_mode,
        latency_s=res.run_start - res.t_send,
        flash_measured_s=[x - res.run_start for x in res.flash_measured],
        rates_before=res.rates_before,
        ambit=dict(name=res.ambit_name, fw=res.ambit_fw),
        felix=res.felix_desc,
        ambit_tags={k: int(len(v)) for k, v in res.ambit.items()},
        saved=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    (out_dir / "run.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log(f"saved to {out_dir}")
    return out_dir


def plot_run(res: RunResult, show: bool = True, save_to: Path | str | None = None):
    import matplotlib.pyplot as plt

    _, edges = time_axis(res.protocol)
    t_amb = res.ambit_time() - res.run_start
    s = res.ambit.get("s_630", np.zeros(0)).astype(float)
    r = res.ambit.get("r_630", np.zeros(0)).astype(float)
    n = min(len(t_amb), len(s), len(r))
    ratio = np.divide(s[:n], r[:n], out=np.full(n, np.nan), where=r[:n] > 0)
    rows = res.summary
    tw = np.array([(x["t_start_s"] + x["t_end_s"]) / 2 for x in rows])
    tot = np.array([x["total"] for x in rows], float)
    tau = np.array([x["tau_mean_ns"] for x in rows], float)
    tau_fit = np.array([x["tau_fit_ns"] for x in rows], float)
    phase = np.array([x["phase"] for x in rows])
    colors = {"pre": "tab:blue", "post": "tab:green", "flash": "tab:red", "edge": "tab:gray", "outside": "lightgray"}

    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=False)
    for ax in axes[:3]:
        for a, b, seg in edges:
            if seg["actinic"] > ACTINIC_ON:
                ax.axvspan(a, b, color="orange", alpha=0.25, lw=0)
    ax = axes[0]
    ax.plot(t_amb[:n], ratio, "o-", ms=3, color="tab:blue")
    ax.set_ylabel("AMBIT s_630 / r_630")
    ax.set_title(
        f"{res.ambit_name}  |  alignment: {res.align_mode}, latency {(res.run_start - res.t_send) * 1e3:.0f} ms"
    )
    ax = axes[1]
    ax.bar(tw, tot, width=res.window_s * 0.9, color=[colors[p] for p in phase])
    ax.set_yscale("log")
    ax.set_ylabel("Felix counts / window")
    ax = axes[2]
    for p in ("pre", "post"):
        m = phase == p
        ax.plot(tw[m], tau[m], "o", color=colors[p], label=f"{p}: mean arrival")
        if np.isfinite(tau_fit[m]).any():
            ax.plot(tw[m], tau_fit[m], "s", mfc="none", color=colors[p], label=f"{p}: tail fit")
    ax.set_ylabel("lifetime proxy [ns]")
    ax.set_xlabel("time since AMBIT sample 0 [s]")
    ax.legend(fontsize=8)
    ax = axes[3]
    t_ns = (res.left_s + (np.arange(res.windows[0].bins.size) + 0.5) * res.bin_s) * 1e9
    for p in ("pre", "post"):
        dec = summed_decay(res, p)
        if dec.sum() > 0:
            ax.semilogy(t_ns, np.maximum(dec, 0.5), color=colors[p], label=f"summed {p} windows")
    ax.set_xlabel("time after SYNC [ns]")
    ax.set_ylabel("counts")
    ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=130)
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--port", help="AMBIT serial port, e.g. COM5 (default: auto-detect)")
    ap.add_argument("--out", help="output folder (default: runs/<timestamp> next to this file)")
    ap.add_argument("--window", type=float, default=WINDOW_S, help="Felix window length [s]")
    ap.add_argument("--threshold", type=float, default=THRESHOLD_V, help="Felix threshold [V]")
    ap.add_argument("--laser-hz", type=float, default=LASER_REP_RATE_HZ)
    ap.add_argument("--pre", type=float, default=PRE_S, help="Felix lead time before arrun [s]")
    ap.add_argument("--post", type=float, default=POST_S, help="Felix tail after the protocol [s]")
    ap.add_argument("--check", action="store_true", help="only detect both instruments")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    ambit = felix = None
    try:
        problems = []
        print("Looking for the AMBIT ...")
        try:
            ambit = find_ambit(args.port)
        except Exception as exc:
            problems.append(f"AMBIT: {exc}")
        print("Looking for the Felix ...")
        try:
            felix = FelixLink().open()
            felix.set_thresholds(args.threshold)
            rates = felix.count_rates(0.5)
            print("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in rates.items()))
            rej = felix.tdc_rejection(0.5)
            print(
                f"Felix TDC drops {rej:.1%} of the SPAD events" if rej == rej else "Felix TDC rejection: not readable"
            )
            warn = rejection_warning(rej, rates.get(SPAD_CH, 0))
            if warn:
                problems.append("Felix: " + warn)
        except Exception as exc:
            problems.append(f"Felix: {exc}")
        if problems:
            for p in problems:
                print("ERROR: " + p)
            return 1
        if args.check:
            return 0
        out = Path(args.out) if args.out else HERE / "runs" / time.strftime("%Y%m%d_%H%M%S")
        res = run_gated(
            ambit,
            felix,
            DEFAULT_PROTOCOL,
            window_s=args.window,
            pre_s=args.pre,
            post_s=args.post,
            laser_hz=args.laser_hz,
            threshold_v=args.threshold,
            out_dir=out,
        )
        plot_run(res, show=not args.no_plot, save_to=out / "summary.png")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1
    finally:
        if felix is not None:
            felix.close()
        if ambit is not None:
            ambit.close()


if __name__ == "__main__":
    sys.exit(main())
