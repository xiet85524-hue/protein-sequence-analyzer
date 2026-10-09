@echo off
set "PY=C:\Users\admin\miniconda3\python.exe"
set "SCRIPT=%~dp0fasta_analysis.py"
set "LOG=%~dp0w1d4_output.log"
echo [BAT-DIR] %~dp0
"%PY%" "%SCRIPT%" all > "%LOG%" 2>&1
echo [EXIT] %errorlevel%
echo [LOG] %LOG%
pause
