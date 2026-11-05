@echo off
rem Launch the phi2 GUI (AMBIT + Felix TCSPC).
rem Double-click it, or from a prompt:  phi2gui.bat --open runs\<run_folder>
rem
rem Runs Python 3.13 (it has PySide6 and pylibtdc), never the conda one: conda's Qt5 DLLs
rem make PySide6 fail with "DLL load failed while importing QtWidgets".
rem Close TDC-Software first: only one program can open the Felix.

setlocal
set "PYTHONPATH="
set "PYTHONHOME="
cd /d "%~dp0scripts\ambit-gate"

set "PY=py -3.13"
where py >nul 2>nul || set "PY=C:\Python313\python.exe"

%PY% -m phi2gui %*
if errorlevel 1 (
    echo.
    echo phi2 GUI stopped with an error, see the messages above.
    pause
)
endlocal
