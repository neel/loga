@echo off
setlocal

rem Directory where this script lives
set "SELF_DIR=%~dp0"
rem Remove trailing backslash
if "%SELF_DIR:~-1%"=="\" set "SELF_DIR=%SELF_DIR:~0,-1%"

rem Prefer DLLs from ./lib for this process
set "PATH=%SELF_DIR%\lib;%PATH%"

rem Run ./bin/loga.exe and forward all args
"%SELF_DIR%\bin\loga.exe" %*

endlocal
