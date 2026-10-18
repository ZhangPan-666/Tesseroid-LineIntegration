@echo off
setlocal
set "REPO=%~dp0.."
cmake -S "%REPO%" -B "%REPO%\build"
if errorlevel 1 exit /b %errorlevel%
cmake --build "%REPO%\build" --config Release --parallel 2
if errorlevel 1 exit /b %errorlevel%
ctest --test-dir "%REPO%\build" -C Release --output-on-failure
exit /b %errorlevel%
