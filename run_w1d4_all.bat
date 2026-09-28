@echo off
chcp 65001 >nul
set "PY=C:\Users\admin\.workbuddy\binaries\python\envs\default\Scripts\python.exe"
set "SCRIPT=C:\Users\admin\Desktop\workbuddy-ai-protein\fasta_analysis.py"
set "LOG=C:\Users\admin\Desktop\workbuddy-ai-protein\w1d4_output.log"
set "CSV=C:\Users\admin\Desktop\workbuddy-ai-protein\w1d4_序列特征对比表.csv"

"%PY%" "%SCRIPT%" all > "%LOG%" 2>&1

if %errorlevel%==0 (
  echo Done. Opening result in Notepad...
) else (
  echo Error %errorlevel%. Open the log to see details.
)
echo Log saved at: %LOG%
echo CSV saved at: %CSV%
echo.
start "" notepad "%LOG%"
echo Press any key to close this window...
pause >nul
