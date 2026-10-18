@echo off
REM Usage: run.bat [processes] [stride] [legacy-test-number]
setlocal
cd /d "%~dp0"
set "EXE=%~dp0..\..\build-mpi\Release\shell_validation.exe"
if not exist "%EXE%" set "EXE=%~dp0..\..\build-mpi\shell_validation.exe"
if not exist "%EXE%" (
    echo Build the shell_validation CMake target first.
    exit /b 1
)
set "NPROC=4"
if not "%~1"=="" set "NPROC=%~1"
set "STRIDE=60"
if not "%~2"=="" set "STRIDE=%~2"
mpiexec -n %NPROC% "%EXE%" %STRIDE% %3
exit /b %errorlevel%
