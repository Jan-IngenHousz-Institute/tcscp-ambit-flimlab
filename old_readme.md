

SPAD detector 9V
Laser 9V
Felix TC 5V

The single-mode fibre P1-S405-FC-1 with small FiberPort PAF2A-2A goes on the laser
The multimode fibre M42L01 with the larger FiberPort PAF2-7B goes on the SPAD


## Sat flash
At 630 nm one photon carries 3.15 × 10⁻¹⁹ J, so one mole of photons is 190 kJ. That gives:

1 µmol s⁻¹ at 630 nm is 0.19 W when this is applied to the surface for the 50 µm fiber core (1.96 × 10⁻⁹ m²) this results in a 3.7 µW on the detector. 

## FLIM GUI (live TCSPC acquisition and lifetime fitting)

Run from the repository root:

    python -m flimgui          # opens the Felix; falls back to the simulator if none is found
    python -m flimgui --sim    # simulator only, no hardware needed

Requirements: Python 3.13 with numpy, scipy, PySide6, pyqtgraph and the pylibtdc wheel
from docs/pylibtdc-2026-09-01.zip (see requirements.txt). Close TDC-Software first: only one
program can open the Felix.

One screen: count-rate tiles for SYNC, CH1 and CH2 (amber above 1 Mcps, red with a beep above
3 Mcps, pile-up ratio shown), threshold and delay per channel (saved on change and on exit,
restored at start), one tab per hardware histogram with multi-hit and accumulate toggles,
IRF recording into irf_library/, and a continuous 1 to 4 exponential fit (reconvolution with the
active IRF, tail fit without one) showing amplitudes, lifetimes, the amplitude-weighted average
lifetime, chi2 and residuals.

Settings live in flimgui_settings.json. Tests: python tests/test_fitting.py


## AMBIT flash protocol with software-gated lifetimes (scripts/ambit-gate)

Runs a PAM-style protocol on an AMBIT leaf sensor (dark, 625 nm actinic flash, relaxation)
while the Felix records contiguous lifetime windows, with no cable between the two:

    cd scripts/ambit-gate
    python ambit_gate.py --check     # detect the AMBIT (serial) and the Felix, print count rates
    python ambit_gate.py             # 1 s dark, 1 s at 250 mA, 3 s relaxation at 10 Hz

The SPAD saturates during the flash, which marks the flash edges in the Felix window totals;
the script aligns the AMBIT samples to those edges and discards the windows that overlap the
flash. The notebook ambit_gate.ipynb walks through the same steps. Why there is no hardware
gate: the AMBIT board has no free output on its connector (GPIO1 drives the flash MOSFET; the
FFC carries only BOOT, Rxd, Txd, reset_EN, GND and V+) and the Felix has no gate input. A
firmware verb pulsing the BOOT line into Felix CH2 would give true per-sample gating; it
needs the Felix input impedance first (50 ohm would hold BOOT low at reset). Details in
scripts/ambit-gate/README.md.
