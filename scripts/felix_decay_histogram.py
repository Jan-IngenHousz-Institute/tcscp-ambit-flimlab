"""Acquire one fluorescence-decay histogram from a TEDIEL Felix TDC with pylibtdc.

Assumed wiring (see docs/ manuals):
    laser SYNC OUT 1  -> Felix SYNC   (TDC channel 0, reference)
    SPAD SMA output   -> Felix CH1    (TDC channel 1, measurement)

Install the API once (Python 3.13 wheel from docs/pylibtdc-2026-09-01.zip):
    pip install --no-index <unzipped>/wheelhouse/pylibtdc-1.18.2.1-cp313-cp313-win_amd64.whl

NOTE: written from the pylibtdc 1.18.2.1 docstrings, not yet run against the
hardware. Confirm the channel numbers with the count-rate printout before
trusting the histogram, and compare with what TDC-Software shows.
"""

import csv
import sys

import pylibtdc as tdc

LASER_REP_RATE_HZ = 20e6          # set on the laser mode selector (5/10/20/40/80 MHz)
WINDOW_S = 1.0 / LASER_REP_RATE_HZ  # 50 ns at 20 MHz: one full laser period
ACQ_SECONDS = 30.0                # integration time of the single histogram
THRESHOLD_V = 1.2                 # both signals are 2.5 V TTL/LVTTL, so ~half
SYNC_CH, SPAD_CH = 0, 1
OUT_CSV = "decay.csv"

RISING = tdc.TDCChannel.edgeTrigger_type.RISING


def main() -> None:
    backend = tdc.TDCBackend()
    devices = backend.getDevices()
    if not devices:
        sys.exit("No Felix found: check the WinUSB driver, USB cable and 5 V supply.")
    print("Opening", devices[0])
    err = backend.open(devices[0])
    print("open() returned:", err)
    backend.loadModules()

    # 1. Input thresholds (leading-edge comparator, 0-2.5 V range).
    thr = backend.getProgrammableThreshold()[0]
    thr.setVoltage([SYNC_CH, SPAD_CH], THRESHOLD_V)

    # 2. Ask the channels for rising edges and keep the requests alive.
    channels = backend.getTDCChannel()
    edge_requests = [channels[c].requestEdge(RISING) for c in (SYNC_CH, SPAD_CH)]

    # 3. Count rates: a sanity check before histogramming.
    counter = backend.getTDCCounter()[0]
    counter.setIntegrationTime(1.0)
    t_int, counts = counter.getCounts()
    # counts[0] is internal; counts[1] is the sync channel; counts[k+1] is channel k.
    sync_rate = counts[SYNC_CH + 1] / t_int
    spad_rate = counts[SPAD_CH + 1] / t_int
    print(f"sync: {sync_rate:,.0f} /s   SPAD: {spad_rate:,.0f} /s")
    if sync_rate < 0.9 * LASER_REP_RATE_HZ:
        print("WARNING: sync rate is below the laser repetition rate. "
              "Check the SYNC cable, the laser mode selector and the threshold.")
    if spad_rate > 0.05 * LASER_REP_RATE_HZ:
        print("WARNING: SPAD rate above 5 % of the laser rate: pile-up. "
              "Lower the laser power or move the collection fibre away.")

    # 4. Hardware histogram: time of each SPAD pulse relative to the last sync.
    hist = backend.getTDCHistogrammer()[0]
    hist.setRefCh(SYNC_CH, RISING)
    hist.setMeasCh(SPAD_CH, RISING)
    hist.setBounds(0.0, WINDOW_S, None)   # seconds; None = smallest bin width
    hist.multiHit = False
    hist.accumulate = True
    left, right, bin_w = hist.getBounds()
    print(f"window {left*1e9:.2f} .. {right*1e9:.2f} ns, bin {bin_w*1e12:.1f} ps")

    print(f"Acquiring for {ACQ_SECONDS:.0f} s ...")
    future = hist.enableSingle(ACQ_SECONDS)
    elapsed, bins = future.get()
    bins = list(bins)
    print(f"done in {elapsed}, total counts {sum(bins):,}, peak bin {max(bins):,}")

    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["time_ns", "counts"])
        for i, c in enumerate(bins):
            writer.writerow([f"{(left + (i + 0.5) * bin_w) * 1e9:.4f}", int(c)])
    print("saved", OUT_CSV)

    del edge_requests
    backend.close()


if __name__ == "__main__":
    main()
