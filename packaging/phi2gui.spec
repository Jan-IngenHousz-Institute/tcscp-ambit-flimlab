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

# One Visual C++ runtime for the whole process. The bundle mixes runtime versions: Python 3.13.5
# brings VCRUNTIME140.dll 14.42 at the top of _internal, PySide6 6.11 brings 14.44 copies of
# MSVCP140*.dll and VCRUNTIME140*.dll in its own folder, and no MSVCP140.dll sits at the top at all.
# Windows loads a DLL once per name: whichever module asks first wins, and a module that finds no
# MSVCP140.dll beside it takes the one of C:\Windows\System32. On a PC whose VC++ redistributable
# is older than Qt 6.11's compiler, Qt then misses a procedure ("DLL load failed while importing
# QtWidgets: The specified procedure could not be found"). Putting PySide6's (newest) copies at the
# top of _internal, which PyInstaller's bootloader puts in the DLL search path before System32,
# makes every module use them, the older importers included (the runtime is backward compatible).
_runtime = ("MSVCP140.dll", "MSVCP140_1.dll", "MSVCP140_2.dll", "VCRUNTIME140.dll", "VCRUNTIME140_1.dll")
_pyside_dir = Path(importlib.metadata.distribution("PySide6").locate_file("PySide6"))
_newest = {name: _pyside_dir / name for name in _runtime if (_pyside_dir / name).is_file()}
if len(_newest) < len(_runtime):
    print(f"WARNING: PySide6 ships only {sorted(_newest)} of the VC runtime; the build keeps Python's copies")
a.binaries = [b for b in a.binaries if not (b[0].lower() in {n.lower() for n in _newest} and "/" not in b[0]
                                               and "\\" not in b[0])]
a.binaries += [(name, str(path), "BINARY") for name, path in _newest.items()]

# Only binaries from this Python and its packages (pylibtdc's loose DLLs live there too).
# PyInstaller resolves the DLLs a binary imports through PATH. Qt6Core.dll of PySide6 6.11 imports
# icuuc.dll, the ICU that Windows 10 (1903+) itself provides in System32 and that PySide6 therefore
# does not ship: built from an Anaconda prompt (Library\bin on PATH), v0.1.3 bundled conda's ICU 73
# instead (icuuc.dll plus a 32 MB icudt73.dll), whose exports carry a version suffix (ucnv_open_73)
# that Qt does not use -> "DLL load failed while importing QtWidgets: The specified procedure could
# not be found" on every PC. The same build took conda's libcrypto/libssl over Python's own.
# So: ICU is never bundled; a binary from outside Python's packages is replaced by Python's own
# copy of that DLL when there is one (DLLs\), else dropped; each case is printed. build_phi2gui.bat
# also clears PATH down to Windows and Python, so normally nothing is found here.
import site
import sys

_roots = [Path(sys.prefix), *map(Path, site.getsitepackages()), Path(site.getusersitepackages()), SRC, ROOT]
_roots = [r.resolve() for r in _roots if r.exists()]
_python_dlls = {p.name.lower(): p for p in [*Path(sys.prefix).glob("*.dll"), *(Path(sys.prefix) / "DLLs").glob("*.dll")]}


def _foreign(src: str) -> bool:
    p = Path(src).resolve()
    return not any(p == r or r in p.parents for r in _roots)


def _icu(dest: str) -> bool:
    return Path(dest).name.lower().startswith(("icuuc", "icuin", "icuio", "icudt", "icutu", "icu."))


_kept = []
for _dest, _src, _kind in a.binaries:
    if _icu(_dest):
        print(f"NOTE: not bundling {_dest} from {_src} (ICU: Windows provides it)")
    elif _foreign(_src):
        _own = _python_dlls.get(Path(_dest).name.lower())
        if _own:
            print(f"NOTE: {_dest}: Python's own {_own} instead of {_src}")
            _kept.append((_dest, str(_own), _kind))
        else:
            print(f"NOTE: not bundling {_dest} from {_src} (outside Python's packages)")
    else:
        _kept.append((_dest, _src, _kind))
a.binaries = _kept

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
