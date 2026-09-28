@echo off
chcp 65001 >nul
set "PY=C:\Users\admin\.workbuddy\binaries\python\versions\3.13.12\python.exe"
set "SCRIPT=C:\Users\admin\Desktop\workbuddy-ai-protein\kgf2_aa.py"
set "LOG=C:\Users\admin\Desktop\workbuddy-ai-protein\w1d1_output.log"

"%PY%" "%SCRIPT%" > "%LOG%" 2>&1

if %errorlevel%==0 (
  echo Done. Opening result in Notepad...
) else (
  echo Error %errorlevel%. Open the log to see details.
)
echo Log saved at: %LOG%
echo.
start "" notepad "%LOG%"
echo Press any key to close this window...
pause >nul
