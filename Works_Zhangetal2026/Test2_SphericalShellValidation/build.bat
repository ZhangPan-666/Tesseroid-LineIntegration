@echo off
setlocal
set "REPO=%~dp0..\.."
cmake -S "%REPO%" -B "%REPO%\build-mpi" -DTLI_BUILD_MPI=ON -DTLI_BUILD_RESEARCH=ON
if errorlevel 1 exit /b %errorlevel%
cmake --build "%REPO%\build-mpi" --config Release --target shell_validation
exit /b %errorlevel%
