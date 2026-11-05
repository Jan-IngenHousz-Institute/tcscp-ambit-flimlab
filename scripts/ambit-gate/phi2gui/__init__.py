"""phi2 GUI: the AMBIT + Felix pauses measurement of notebook 4 with its results, and the upload
of a run to an openJII experiment.

Run it from anywhere::

    python -m phi2gui                 # from scripts/ambit-gate
    python -m phi2gui --open <run>    # show a saved run folder

The package lives next to the ambit-gate modules (ambit_gate, ambit_pause, reconv, ...) and puts
that folder on sys.path, so it needs no installation.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))
