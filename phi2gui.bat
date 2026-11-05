@echo off
rem Launch the phi2 GUI (AMBIT + Felix TCSPC) from a source checkout of the repository.
rem Double-click it, or from a prompt:  phi2gui.bat --open runs\<run_folder>
rem
rem Needs Python 3.13 with the packages of requirements.txt and the pylibtdc wheel installed.
rem No Python on this PC? Use the standalone build instead: phi2gui-v<version>-win64.zip
rem (on the GitHub Releases page, built with build_phi2gui.bat); it needs nothing installed.
rem
rem Runs Python 3.13 (it has PySide6 and pylibtdc), never the conda one: conda's Qt5 DLLs
rem make PySide6 fail with "DLL load failed while importing QtWidgets".
rem Close TDC-Software first: only one program can open the Felix.

setlocal
set "PYTHONPATH="
set "PYTHONHOME="

if not exist "%~dp0scripts\ambit-gate\phi2gui\__main__.py" (
    echo phi2gui.bat must stay in the repository folder, next to scripts\ambit-gate.
    echo If you opened it from inside a zip file, extract the whole zip first.
    echo Looked for: %~dp0scripts\ambit-gate\phi2gui
    goto :fail
)
cd /d "%~dp0scripts\ambit-gate"

set "PY="
where py >nul 2>nul && py -3.13 -c "pass" >nul 2>nul && set "PY=py -3.13"
if not defined PY if exist "C:\Python313\python.exe" set "PY=C:\Python313\python.exe"
if not defined PY (
    echo Python 3.13 was not found on this PC.
    echo Either install it from https://www.python.org/downloads/ and run
    echo     pip install -r requirements.txt   ^(plus the pylibtdc wheel, see README.md^)
    echo or use the standalone build phi2gui-v^<version^>-win64.zip, which needs no Python.
    goto :fail
)

%PY% -m phi2gui %*
if errorlevel 1 goto :fail
endlocal
exit /b 0

:fail
echo.
echo phi2 GUI stopped with an error, see the messages above.
pause
endlocal
exit /b 1
