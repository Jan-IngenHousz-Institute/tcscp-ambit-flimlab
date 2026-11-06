"""Felix WinUSB driver: find the bundled TEDIEL installer and run it with administrator rights.

``Felix_WinUSB_driver.exe`` is TEDIEL's free download (https://tediel.com/resources), the same file
TDC-Software's setup runs. It extracts the driver package, signs its catalogue (no "unverified
publisher" warning) and runs Microsoft's dpinst: Windows asks for administrator permission once,
then a two-click wizard installs the WinUSB driver for the Felix (USB VID 0403, PID 6C68).
"""

from __future__ import annotations

import sys
from pathlib import Path

INSTALLER_NAME = "Felix_WinUSB_driver.exe"
DOWNLOAD_URL = "https://tediel.com/wp-content/Resources/Felix/Felix_WinUSB_driver.exe"
ERROR_CANCELLED = 1223  # the UAC prompt was refused


def installer_path() -> Path | None:
    """The bundled installer: ``_internal/driver/`` in the frozen app, ``docs/`` in a source checkout."""
    if getattr(sys, "frozen", False):
        candidates = [
            Path(getattr(sys, "_MEIPASS", "")) / "driver" / INSTALLER_NAME,
            Path(sys.executable).resolve().parent / INSTALLER_NAME,
        ]
    else:
        candidates = [Path(__file__).resolve().parents[3] / "docs" / INSTALLER_NAME]
    return next((p for p in candidates if p.is_file()), None)


def run_installer(path: Path) -> int:
    """Run the installer elevated (UAC prompt) and wait for it; returns its exit code.

    Raises ``PermissionError`` when the UAC prompt is refused, ``OSError`` when it cannot start.
    """
    if sys.platform != "win32":
        raise OSError("the Felix driver installer is for Windows only")
    import ctypes
    import ctypes.wintypes as wt

    class ShellExecuteInfo(ctypes.Structure):
        _fields_ = [
            ("cbSize", wt.DWORD),
            ("fMask", wt.ULONG),
            ("hwnd", wt.HWND),
            ("lpVerb", wt.LPCWSTR),
            ("lpFile", wt.LPCWSTR),
            ("lpParameters", wt.LPCWSTR),
            ("lpDirectory", wt.LPCWSTR),
            ("nShow", ctypes.c_int),
            ("hInstApp", wt.HINSTANCE),
            ("lpIDList", ctypes.c_void_p),
            ("lpClass", wt.LPCWSTR),
            ("hkeyClass", wt.HKEY),
            ("dwHotKey", wt.DWORD),
            ("hIcon", wt.HANDLE),
            ("hProcess", wt.HANDLE),
        ]

    shell32 = ctypes.WinDLL("shell32", use_last_error=True)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.WaitForSingleObject.argtypes = [wt.HANDLE, wt.DWORD]
    kernel32.GetExitCodeProcess.argtypes = [wt.HANDLE, ctypes.POINTER(wt.DWORD)]
    kernel32.CloseHandle.argtypes = [wt.HANDLE]
    info = ShellExecuteInfo(
        cbSize=ctypes.sizeof(ShellExecuteInfo),
        fMask=0x40,  # SEE_MASK_NOCLOSEPROCESS
        lpVerb="runas",
        lpFile=str(path),
        lpParameters="",
        lpDirectory=str(path.parent),
        nShow=1,
    )
    if not shell32.ShellExecuteExW(ctypes.byref(info)):
        err = ctypes.get_last_error()
        if err == ERROR_CANCELLED:
            raise PermissionError("administrator permission was refused, the driver was not installed")
        raise OSError(err, f"could not start {path.name} (Windows error {err})")
    if not info.hProcess:
        return 0
    kernel32.WaitForSingleObject(info.hProcess, 0xFFFFFFFF)
    code = wt.DWORD()
    kernel32.GetExitCodeProcess(info.hProcess, ctypes.byref(code))
    kernel32.CloseHandle(info.hProcess)
    return int(code.value)
