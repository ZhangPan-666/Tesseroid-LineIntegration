@echo off
setlocal
set "REPO=%~dp0.."
cmake -S "%REPO%" -B "%REPO%\build-research" -DTLI_BUILD_RESEARCH=ON
if errorlevel 1 exit /b %errorlevel%
cmake --build "%REPO%\build-research" --config Release --parallel 2 --target kernel_peak_regimes peak_location_study radial_peak_profiles
exit /b %errorlevel%
