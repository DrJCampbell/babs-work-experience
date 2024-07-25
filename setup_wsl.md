# Setting up WSL for work experience laptops

We are unable to get loan Mac laptops for work experience students and so are using Dell laptops running Microsoft Windows.

Windows enables installation of a Linux distrubution through the Windows Subsystem for Linux (WSL). These steps were followed into set up WSL with Ubuntu. This should work with Windows 10 but was found to be unreliable. The steps below worked on Windows 11.

1. Click the Start button on the Windows Task Bar.
2. Find the CMD.exe application, right-click it and choose 'Run as Administrator'
3. In the CMD.exe window, run the command `wsl --install`
4. If you hit a problem, try the following:

```sh
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

All being well, you should see an orange Ubuntu icon appear on the start menu.


