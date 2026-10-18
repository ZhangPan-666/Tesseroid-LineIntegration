@echo off
setlocal

if not defined TLI_ONEAPI_ROOT set "TLI_ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI"
if not defined TLI_MATLAB_ROOT (
    echo Set TLI_MATLAB_ROOT to your MATLAB installation directory first.
    exit /b 1
)
set "FUNCTIONS_DIR=%~dp0"
for %%I in ("%FUNCTIONS_DIR%..\C++Codes") do set "CPP_DIR=%%~fI"

call "%TLI_ONEAPI_ROOT%\setvars.bat" >nul 2>&1
if errorlevel 1 (
    echo Failed to initialize Intel oneAPI.
    exit /b 1
)

icx-cl /nologo /O2 /std:c++17 /EHsc /MD /Qopenmp ^
    /DMATLAB_MEX_FILE /DMW_NEEDS_VERSION_H ^
    /I"%TLI_MATLAB_ROOT%\extern\include" /I"%CPP_DIR%" ^
    /c "%FUNCTIONS_DIR%tli_lineintegration_mex.cpp" ^
    /Fo:"%FUNCTIONS_DIR%tli_lineintegration_mex.obj"
if errorlevel 1 (
    echo MEX compile failed.
    exit /b 1
)

icx-cl /nologo /LD /MD /Qopenmp ^
    "%FUNCTIONS_DIR%tli_lineintegration_mex.obj" ^
    /link /OUT:"%FUNCTIONS_DIR%tli_lineintegration_mex.mexw64" ^
    /LIBPATH:"%TLI_MATLAB_ROOT%\extern\lib\win64\microsoft" ^
    libmex.lib libmx.lib /EXPORT:mexFunction

if errorlevel 1 (
    echo MEX build failed.
    exit /b 1
)

del /q "%FUNCTIONS_DIR%tli_lineintegration_mex.exp" 2>nul
del /q "%FUNCTIONS_DIR%tli_lineintegration_mex.lib" 2>nul
del /q "%FUNCTIONS_DIR%tli_lineintegration_mex.obj" 2>nul

echo Built "%FUNCTIONS_DIR%tli_lineintegration_mex.mexw64"
exit /b 0
