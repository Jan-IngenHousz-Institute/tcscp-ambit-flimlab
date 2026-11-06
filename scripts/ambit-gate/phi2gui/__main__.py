"""Entry point: python -m phi2gui [--open RUN_FOLDER] [--settings FILE] [--smoke-test FILE]"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path


def smoke_test(out: Path) -> int:
    """Load everything a real start loads (Qt, pylibtdc and its DLLs, the main window) without
    showing anything, then write ``out``. The release workflow runs it on the frozen build, where
    a missing or foreign DLL ("DLL load failed while importing QtWidgets") stops the release."""
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    if sys.platform != "darwin":  # TEDIEL ships pylibtdc for Windows and Linux only
        import pylibtdc  # noqa: F401  (the Felix bindings load lazily in a real run)
    from PySide6.QtWidgets import QApplication

    from . import __version__
    from . import openjii as oj
    from .app import MainWindow

    app = QApplication(sys.argv[:1])
    with tempfile.TemporaryDirectory() as tmp:
        win = MainWindow(oj.Settings(Path(tmp) / "phi2gui.json"))
        app.processEvents()
        win.deleteLater()
    app.processEvents()
    out.write_text(f"ok {__version__}\n", encoding="utf-8")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="AMBIT + Felix phi2 measurement with openJII upload")
    parser.add_argument("--open", help="show this saved run folder at start (no instruments needed)")
    parser.add_argument("--settings", help="settings file (default: the user's profile)")
    parser.add_argument("--smoke-test", metavar="FILE", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.smoke_test:
        return smoke_test(Path(args.smoke_test))

    from PySide6.QtWidgets import QApplication

    from . import openjii as oj
    from .app import MainWindow

    app = QApplication(sys.argv)
    app.setApplicationName("phi2 GUI")
    win = MainWindow(oj.Settings(Path(args.settings) if args.settings else None))
    win.show()
    if args.open:
        win.open_run(Path(args.open))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
