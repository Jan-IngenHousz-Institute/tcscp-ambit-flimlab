# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the standalone phi2 GUI.

Build with build_phi2gui.bat in the repo root (or: py -3.13 -m PyInstaller packaging/phi2gui.spec).
Result: dist/phi2gui/phi2gui.exe plus an _internal folder; zip the phi2gui folder to share it.

The bundle includes TEDIEL's pylibtdc (Felix TDC bindings, authored by Politecnico di Milano),
the same wheel TEDIEL offers as a free download at https://tediel.com/resources.
"""
import importlib.metadata
from pathlib import Path

ROOT = Path(SPECPATH).resolve().parent          # repo root (this file is in packaging/)
SRC = ROOT / "scripts" / "ambit-gate"           # phi2gui and the modules it imports as siblings

# pylibtdc is one .pyd whose DLLs (TDC, libusb, ftdi1, ...) are dropped loose into site-packages.
# Take every DLL the wheel installed so the Felix works on a PC without Python.
_tdc = importlib.metadata.distribution("pylibtdc")
pylibtdc_dlls = [(str(_tdc.locate_file(f)), ".") for f in _tdc.files if f.suffix.lower() == ".dll"]

# TEDIEL's Felix WinUSB driver installer: the GUI offers to run it when the Felix has no driver yet
# (build_phi2gui.bat downloads it into docs/ when missing).
_drv = ROOT / "docs" / "Felix_WinUSB_driver.exe"
if not _drv.is_file():
    print(f"WARNING: {_drv} missing: this build cannot install the Felix driver")
datas = [(str(_drv), "driver")] if _drv.is_file() else []
# the IRF library (blank_dose.DEFAULT_IRF_NAME is the GUI's default; users add their own via the GUI)
datas += [(str(p), "irf-library") for p in sorted((ROOT / "irf-library").glob("*.csv"))]

a = Analysis(
    [str(ROOT / "packaging" / "phi2gui_entry.py")],
    pathex=[str(SRC), str(ROOT)],
    binaries=pylibtdc_dlls,
    datas=datas,
    hiddenimports=[
        "pylibtdc",                         # imported lazily inside ambit_gate (open_felix)
        "flimgui.fitting",                  # optional tail fit, imported from the repo root
        "serial.tools.list_ports",
        "matplotlib.backends.backend_agg",
        "matplotlib.backends.backend_qtagg",
    ],
    excludes=[
        "tkinter", "IPython", "jupyter", "notebook", "pytest",
        "pandas", "sklearn", "xgboost", "statsmodels", "seaborn", "h5py", "pyarrow",
        "pyqtgraph", "PyQt5", "PyQt6",
        # pulled in through optional imports of requests/urllib3, never loaded by phi2gui (~20 MB)
        "cryptography", "OpenSSL", "paramiko", "bcrypt", "nacl", "jinja2", "psutil",
    ],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    exclude_binaries=True,
    name="phi2gui",
    console=False,                          # GUI app: no console window; stdout/stderr go to phi2gui.log
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name="phi2gui", upx=False)
