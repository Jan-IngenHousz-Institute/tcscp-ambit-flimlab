"""PyInstaller entry point for the phi2 GUI (see packaging/phi2gui.spec and build_phi2gui.bat).

The frozen build is a windowed app without a console, so whatever the program writes to
stdout/stderr (tracebacks of background tasks, Qt warnings) is kept in phi2gui.log next to
phi2gui.exe. Running from source (phi2gui.bat / python -m phi2gui) is unchanged.
"""
import sys
from pathlib import Path

if getattr(sys, "frozen", False):
    _log = open(Path(sys.executable).with_name("phi2gui.log"), "a", encoding="utf-8", buffering=1)
    sys.stdout = sys.stderr = _log

from phi2gui.__main__ import main  # noqa: E402

sys.exit(main())
