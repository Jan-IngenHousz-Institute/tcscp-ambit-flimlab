"""PyInstaller entry point for the phi2 GUI (see packaging/phi2gui.spec and build_phi2gui.bat).

The frozen build is a windowed app without a console, so whatever the program writes to
stdout/stderr (tracebacks of background tasks, Qt warnings) is kept in phi2gui.log: next to the
executable on Windows and Linux, in ~/Library/Logs on macOS, where the .app bundle must stay
unmodified. Running from source (phi2gui.bat / python -m phi2gui) is unchanged.
"""

import sys
from pathlib import Path


def log_path() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Logs" / "phi2gui.log"
    return Path(sys.executable).with_name("phi2gui.log")


if getattr(sys, "frozen", False):
    _path = log_path()
    _path.parent.mkdir(parents=True, exist_ok=True)
    _log = open(_path, "a", encoding="utf-8", buffering=1)
    sys.stdout = sys.stderr = _log

from phi2gui.__main__ import main  # noqa: E402

sys.exit(main())
