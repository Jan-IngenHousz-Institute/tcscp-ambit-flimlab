# Documentation

The repository's [README](../README.md) is for using the set-up: wiring, the measurement, the phi2
GUI, installing. The pages here are for changing it.

| Page | Read it when |
| --- | --- |
| [architecture.md](architecture.md) | You change the library, the phi2 GUI or flimgui, and need to know how the pieces fit, what a run folder holds, or where settings live |
| [openjii-upload.md](openjii-upload.md) | You touch what the phi2 GUI uploads, or test an upload against openJII |
| [releasing.md](releasing.md) | You want a new standalone build, change the packaging, or move to a new pylibtdc |

Agents start at [AGENTS.md](../AGENTS.md).

## Vendor manuals

| File | What it is |
| --- | --- |
| `TDC-SW_manual.pdf` | TEDIEL's TDC Software user guide for the Felix |
| `Histogram_FastTimeTrigger.pdf` | The Felix hardware histogram (what `pylibtdc`'s histogram module computes) |
| `binaryFormatGuide.pdf` | libTDC's TTM binary format (time tags), which `ambit_chop` and `chop_sim` read and write |
| `PLD-manual.pdf` | FLIM LABS picosecond pulsed fibre-coupled laser |
| `SPAD-manual.pdf` | FLIM LABS single-photon SPAD detector SPAD-050-C-B-FC |
| `ambit-tcscp-adapter.png` | The AMBIT-to-TCSPC adapter, shown in the README |

TEDIEL's installers are not committed. `packaging/fetch_pylibtdc.py` puts
`pylibtdc-2026-09-01.zip` here, and `build_phi2gui.bat` puts `Felix_WinUSB_driver.exe` here, both
from https://tediel.com/resources.
