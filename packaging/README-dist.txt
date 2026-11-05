phi2 GUI  -  AMBIT + Felix TCSPC phi2 measurement with openJII upload
=====================================================================

Run:      double-click phi2gui.exe  (no Python needed)
          phi2gui.exe --open <run folder>   shows a saved run, no instruments needed
          If Windows SmartScreen warns, click "More info" > "Run anyway" (build is not code-signed).

Before the first measurement on a new PC:
  1. Plug in the Felix and the AMBIT (USB-serial, 115200 8N1), then start phi2gui.exe.
  2. If Windows has no driver for the Felix yet, the app offers to install TEDIEL's WinUSB
     driver: allow the administrator prompt, then Next and Finish in the short wizard.
     (Manual alternative: Felix_WinUSB_driver.exe from https://tediel.com/resources.)
  3. Close TDC-Software: only one program can open the Felix.

IRF:      the instrument response for the lifetime fits is shipped (261104 black plastic, 50 %);
          "Add IRF file..." takes a TDC-Software histogram export (CSV) and remembers it.
Rule:     the Felix window must be at most half of the shortest pause (the app checks before a run).
Output:   runs\<timestamp>_gui\ next to phi2gui.exe (raw.png shows the raw data of every run)
Settings: %APPDATA%\phi2gui\phi2gui.json  (openJII API key, last used ports)
Log:      phi2gui.log next to phi2gui.exe  (send it along when reporting a problem)

Keep this folder together: phi2gui.exe needs the _internal folder beside it.
This package bundles TEDIEL's pylibtdc library for the Felix TDC (free download at
https://tediel.com/resources, Python API 2026.09.01); it only works with a licensed Felix.

Source and releases: https://github.com/Jan-IngenHousz-Institute/tcscp-ambit-flimlab
