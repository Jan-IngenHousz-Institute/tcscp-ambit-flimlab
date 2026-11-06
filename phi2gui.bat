@echo off
rem Launch the phi2 GUI (AMBIT + Felix TCSPC) from a source checkout of the repository.
rem Double-click it, or from a prompt:  phi2gui.bat --open runs\<run_folder>
rem
rem With uv installed (winget install astral-sh.uv) it sets everything up itself the first time:
rem Python 3.13, the packages of uv.lock and TEDIEL's pylibtdc, in .venv next to this file.
rem Without uv it falls back to an installed Python 3.13 that already has the packages.
rem No Python on this PC? Use the standalone build instead: phi2gui-v<version>-win64.zip
rem on the GitHub Releases page; it needs nothing installed.
rem
rem Never the conda Python: conda's Qt5 DLLs make PySide6 fail with
rem "DLL load failed while importing QtWidgets". Close TDC-Software first: only one program can
rem open the Felix.

setlocal
set "PYTHONPATH="
set "PYTHONHOME="

if not exist "%~dp0scripts\ambit-gate\phi2gui\__main__.py" (
    echo phi2gui.bat must stay in the repository folder, next to scripts\ambit-gate.
    echo If you opened it from inside a zip file, extract the whole zip first.
    echo Looked for: %~dp0scripts\ambit-gate\phi2gui
    goto :fail
)

set "UV="
for /f "delims=" %%i in ('where uv 2^>nul') do if not defined UV set "UV=%%i"
if defined UV (
    pushd "%~dp0"
    "%UV%" sync --locked --quiet || (popd & goto :fail)
    popd
    cd /d "%~dp0scripts\ambit-gate"
    "%~dp0.venv\Scripts\python.exe" -m phi2gui %*
    if errorlevel 1 goto :fail
    endlocal
    exit /b 0
)

cd /d "%~dp0scripts\ambit-gate"
set "PY="
where py >nul 2>nul && py -3.13 -c "pass" >nul 2>nul && set "PY=py -3.13"
if not defined PY if exist "C:\Python313\python.exe" set "PY=C:\Python313\python.exe"
if not defined PY (
    echo Neither uv nor Python 3.13 was found on this PC.
    echo Install uv with  winget install astral-sh.uv  and double-click phi2gui.bat again,
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
