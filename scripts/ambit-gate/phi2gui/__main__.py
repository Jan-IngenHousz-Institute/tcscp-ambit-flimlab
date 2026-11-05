"""Entry point: python -m phi2gui [--open RUN_FOLDER] [--settings FILE]"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="AMBIT + Felix phi2 measurement with openJII upload")
    parser.add_argument("--open", help="show this saved run folder at start (no instruments needed)")
    parser.add_argument("--settings", help="settings file (default: the user's profile)")
    args = parser.parse_args(argv)

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
