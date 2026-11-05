phi2 GUI  -  AMBIT + Felix TCSPC phi2 measurement with openJII upload
=====================================================================

Run:      double-click phi2gui.exe  (no Python needed)
          phi2gui.exe --open <run folder>   shows a saved run, no instruments needed

Before the first measurement on a new PC:
  1. Install the Felix WinUSB driver (Felix_WinUSB_driver.exe from FLIM LABS), once.
  2. Plug in the Felix and the AMBIT (USB-serial, 115200 8N1).
  3. Close TDC-Software: only one program can open the Felix.

Output:   runs\<timestamp>_gui\ next to phi2gui.exe
Settings: %APPDATA%\phi2gui\phi2gui.json  (openJII API key, last used ports)
Log:      phi2gui.log next to phi2gui.exe  (send it along when reporting a problem)

Keep this folder together: phi2gui.exe needs the _internal folder beside it.
This package contains the FLIM LABS pylibtdc library; for internal use only.

Source: https://github.com/Jan-IngenHousz-Institute/tcscp-ambit-flimlab
