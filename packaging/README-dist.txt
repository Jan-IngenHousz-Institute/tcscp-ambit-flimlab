phi2 GUI: AMBIT + Felix TCSPC phi2 measurement with openJII upload
===================================================================

This folder holds everything the app needs. You install nothing else, not even Python.
Keep the folder together: the program needs the files beside it.

Windows 10 or 11, 64-bit
  Start:     double-click phi2gui.exe
  Warning:   the build is not code-signed. If SmartScreen warns, click "More info", then "Run anyway".
  Felix:     if Windows has no driver for the Felix yet, the app offers to install TEDIEL's WinUSB
             driver. Allow the administrator prompt, then click Next and Finish in the short wizard.
             Felix_WinUSB_driver.exe from https://tediel.com/resources does the same by hand.
  Output:    runs\<timestamp>_gui\ next to phi2gui.exe
  Log:       phi2gui.log next to phi2gui.exe
  Settings:  %APPDATA%\phi2gui\phi2gui.json, with your openJII API key and the last used ports

Linux x64, any distribution with glibc 2.35 or newer
  Start:     ./phi2gui in this folder
  Felix:     a normal user needs permission for the USB device. As root, put this line in
             /etc/udev/rules.d/70-felix.rules, then unplug and replug the Felix:
               SUBSYSTEM=="usb", ATTRS{idVendor}=="0403", ATTRS{idProduct}=="6c68", TAG+="uaccess"
             For the AMBIT's serial port, add yourself to the dialout group.
             Nobody has measured with the Linux build yet. Report what you find.
  Output:    runs/<timestamp>_gui/ next to phi2gui
  Log:       phi2gui.log next to phi2gui
  Settings:  ~/.config/phi2gui/phi2gui.json

macOS, Apple silicon
  Start:     move phi2gui.app to Applications and double-click it. The app is not notarized, so
             macOS blocks the first start. Click Done, open System Settings > Privacy & Security,
             scroll to Security, click "Open Anyway" next to phi2gui and confirm with your password.
             From then on it starts with a double-click.
  Damaged:   if macOS instead says the app is damaged and offers to move it to the Bin, run this
             once in Terminal, then start it as above:
               xattr -dr com.apple.quarantine /Applications/phi2gui.app
  Felix:     TEDIEL ships no macOS library, so a Mac cannot drive the Felix. It opens saved runs,
             shows their results and uploads them to openJII.
  Output:    ~/Documents/phi2gui/runs/
  Log:       ~/Library/Logs/phi2gui.log
  Settings:  ~/Library/Application Support/phi2gui/phi2gui.json

Every computer
  "Open saved run..." shows a run folder again without instruments, on any of the three.
  Close TDC-Software before measuring. Only one program can open the Felix at a time.
  The lifetime fits use the shipped instrument response, 261104 black plastic at 50 %. "Add IRF
  file..." takes a TDC-Software histogram export in CSV and remembers it.
  The Felix window must be at most half of the shortest pause. The app checks this before a run.
  Each run folder has raw.png with the raw data, also when the fits fail.
  Send phi2gui.log along when you report a problem.

The app bundles TEDIEL's pylibtdc library, Python API 2026.09.01, which is a free download at
https://tediel.com/resources. It only works with a licensed Felix.

Source and releases: https://github.com/Jan-IngenHousz-Institute/tcscp-ambit-flimlab
