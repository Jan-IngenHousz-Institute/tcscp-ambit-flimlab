# Documentation

The repository's [README](../README.md) covers using the set-up: wiring, the measurement, the phi2
GUI and installing it. The pages here cover changing it.

| Page | Read it when |
| --- | --- |
| [architecture.md](architecture.md) | You change the library, the phi2 GUI or flimgui and need to know how the pieces fit, what a run folder holds or where settings live |
| [openjii-upload.md](openjii-upload.md) | You touch what the phi2 GUI uploads, or test an upload against openJII |
| [releasing.md](releasing.md) | You want a new standalone build, change the packaging, or move to a new pylibtdc |

Agents start at [AGENTS.md](../AGENTS.md).

## Vendor manuals

| File | What it is |
| --- | --- |
| `TDC-SW_manual.pdf` | TEDIEL's user guide for TDC Software, which drives the Felix |
| `Histogram_FastTimeTrigger.pdf` | The Felix hardware histogram, which `pylibtdc`'s histogram module computes |
| `binaryFormatGuide.pdf` | libTDC's TTM binary format for time tags, which `ambit_chop` reads and `chop_sim` writes |
| `PLD-manual.pdf` | The FLIM LABS picosecond pulsed fibre-coupled laser |
| `SPAD-manual.pdf` | The FLIM LABS SPAD-050-C-B-FC single-photon detector |
| `ambit-tcscp-adapter.png` | The AMBIT-to-TCSPC adapter shown in the README |

TEDIEL's installers stay out of git. `build_phi2gui.bat` and the release workflow download
`Felix_WinUSB_driver.exe` into this folder from https://tediel.com/resources.
