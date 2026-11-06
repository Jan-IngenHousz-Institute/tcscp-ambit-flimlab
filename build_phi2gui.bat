@echo off
rem Build the standalone phi2 GUI locally and zip it:
rem     dist\phi2gui\phi2gui.exe            (plus _internal\, README.txt)
rem     dist\phi2gui-v<version>-win64.zip   (version = phi2gui.__version__)
rem
rem Releases do not need this: .github\workflows\gui-release.yml builds and publishes one on every
rem merge to main. This is for trying a build on your own PC.
rem
rem Needs uv (https://docs.astral.sh/uv/, e.g.  winget install astral-sh.uv). uv brings Python 3.13
rem and every package from uv.lock into .venv; packaging\fetch_pylibtdc.py puts TEDIEL's pylibtdc
rem (free download at https://tediel.com/resources) into wheelhouse\ first.

setlocal
set "PYTHONPATH="
set "PYTHONHOME="
cd /d "%~dp0"

set "UV="
for /f "delims=" %%i in ('where uv 2^>nul') do if not defined UV set "UV=%%i"
if not defined UV (
    echo uv was not found. Install it with:  winget install astral-sh.uv
    goto :fail
)

"%UV%" run --no-project packaging\fetch_pylibtdc.py || goto :fail
"%UV%" sync --locked --group build || goto :fail
set "VPY=%~dp0.venv\Scripts\python.exe"
"%VPY%" -c "import sys; print(sys.base_prefix)" > "%TEMP%\phi2gui_base.txt" || goto :fail
set /p BASEPY=<"%TEMP%\phi2gui_base.txt"

"%VPY%" -c "import sys; sys.path.insert(0, r'scripts\ambit-gate'); import phi2gui; print(phi2gui.__version__)" > "%TEMP%\phi2gui_version.txt" || goto :fail
set /p VER=<"%TEMP%\phi2gui_version.txt"
if "%VER%"=="" goto :fail

if not exist docs\Felix_WinUSB_driver.exe (
    echo Fetching TEDIEL's Felix WinUSB driver installer into docs\ ...
    curl -L -o docs\Felix_WinUSB_driver.exe https://tediel.com/wp-content/Resources/Felix/Felix_WinUSB_driver.exe || goto :fail
)

rem A PATH of Windows and this Python only. PyInstaller resolves the DLLs a binary imports through
rem PATH, and from an Anaconda prompt it bundled conda's ICU over the one Windows provides for Qt
rem (v0.1.3 did not start on any PC). packaging\phi2gui.spec drops such foreign DLLs as well.
set "PATH=%SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\WindowsPowerShell\v1.0;%BASEPY%;%~dp0.venv\Scripts"

echo Building phi2gui v%VER% ...
"%VPY%" -m PyInstaller --noconfirm --clean --distpath dist --workpath build packaging\phi2gui.spec || goto :fail
copy /y packaging\README-dist.txt dist\phi2gui\README.txt >nul

rem The exe loads Qt, pylibtdc's DLLs and the main window, then exits (as the release workflow checks).
if exist "%TEMP%\phi2gui_smoke.txt" del "%TEMP%\phi2gui_smoke.txt"
start "" /wait dist\phi2gui\phi2gui.exe --smoke-test "%TEMP%\phi2gui_smoke.txt"
findstr /b /c:"ok %VER%" "%TEMP%\phi2gui_smoke.txt" >nul 2>nul || (
    echo The built exe did not start, see dist\phi2gui\phi2gui.log
    goto :fail
)
del dist\phi2gui\phi2gui.log 2>nul

set "ZIP=dist\phi2gui-v%VER%-win64.zip"
if exist "%ZIP%" del "%ZIP%"
powershell -NoProfile -Command "Compress-Archive -Path 'dist\phi2gui' -DestinationPath '%ZIP%' -CompressionLevel Optimal" || goto :fail

echo.
echo Built %ZIP%
endlocal
exit /b 0

:fail
echo.
echo Build failed, see the messages above.
endlocal
exit /b 1
