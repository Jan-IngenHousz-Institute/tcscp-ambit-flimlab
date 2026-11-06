"""Entry point: python -m flimgui [--sim]"""

from __future__ import annotations

import argparse
import logging
import sys

from PySide6.QtWidgets import QApplication, QMessageBox

from .device import DeviceError, FelixDevice, SimDevice
from .irf import IRFStore
from .main_window import MainWindow
from .settings import Settings


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Live TCSPC acquisition and lifetime fitting (Felix TDC)")
    parser.add_argument("--sim", action="store_true", help="run with the simulator, no hardware")
    parser.add_argument("--serial", default="", help="serial number of the Felix to open")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    app = QApplication(sys.argv)
    app.setApplicationName("FLIM GUI")
    settings = Settings.load()

    device = None
    if not args.sim:
        device = FelixDevice(args.serial)
        try:
            device.open()
        except DeviceError as exc:
            QMessageBox.warning(None, "No Felix", f"{exc}\n\nStarting in simulator mode so you can still use the GUI.")
            device = None
        except Exception as exc:  # unexpected library error
            QMessageBox.warning(
                None, "Device error", f"Unexpected error while opening the Felix:\n{exc}\n\nStarting in simulator mode."
            )
            device = None
    if device is None:
        device = SimDevice()
        device.open()
    settings.sim_mode = device.is_simulator

    win = MainWindow(device, settings, IRFStore())
    win.resize(1500, 900)
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
