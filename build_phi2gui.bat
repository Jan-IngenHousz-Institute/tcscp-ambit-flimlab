@echo off
rem Build the standalone phi2 GUI and zip it for colleagues:
rem     dist\phi2gui\phi2gui.exe            (plus _internal\, README.txt)
rem     dist\phi2gui-v<version>-win64.zip   (version = phi2gui.__version__)
rem
rem Needs Python 3.13 with requirements.txt, pylibtdc (wheel in docs\pylibtdc-*.zip) and
rem pyinstaller installed. Same Python rule as phi2gui.bat: never the conda one.
rem
rem The zip contains the FLIM LABS pylibtdc library: share it internally (OneDrive, Teams),
rem do not attach it to the public GitHub repo.

setlocal
set "PYTHONPATH="
set "PYTHONHOME="
cd /d "%~dp0"

set "PY=py -3.13"
where py >nul 2>nul || set "PY=C:\Python313\python.exe"

%PY% -c "import sys; sys.path.insert(0, r'scripts\ambit-gate'); import phi2gui; print(phi2gui.__version__)" > "%TEMP%\phi2gui_version.txt" || goto :fail
set /p VER=<"%TEMP%\phi2gui_version.txt"
if "%VER%"=="" goto :fail

echo Building phi2gui v%VER% ...
%PY% -m PyInstaller --noconfirm --clean --distpath dist --workpath build packaging\phi2gui.spec || goto :fail

copy /y packaging\README-dist.txt dist\phi2gui\README.txt >nul

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
