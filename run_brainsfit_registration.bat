@echo off
setlocal

REM ============================================================
REM BRAINSFit batch registration launcher
REM
REM Instructions:
REM 1. Set SLICER_EXE to the location of your 3D Slicer executable.
REM 2. Keep this file in the same folder as brainsfit_batch_registration.py.
REM 3. Run this batch file to start the registration.
REM
REM The existing registration results will not be overwritten.
REM ============================================================

set "SLICER_EXE=C:\Path\To\3D Slicer\Slicer.exe"
set "REGISTRATION_SCRIPT=%~dp0brainsfit_batch_registration.py"

if not exist "%SLICER_EXE%" (
    echo ERROR: 3D Slicer was not found at:
    echo %SLICER_EXE%
    echo.
    echo Open 3D Slicer, select View - Python Interactor, and run:
    echo print^(slicer.app.applicationFilePath^(^)^)
    echo Then replace SLICER_EXE in this file with the printed path.
    pause
    exit /b 1
)

if not exist "%REGISTRATION_SCRIPT%" (
    echo ERROR: Registration script was not found at:
    echo %REGISTRATION_SCRIPT%
    pause
    exit /b 1
)

echo Starting BRAINSFit batch registration.
echo The existing registration results will not be overwritten.
echo This may take a long time for all datasets.
echo.

"%SLICER_EXE%" --no-splash --no-main-window --python-script "%REGISTRATION_SCRIPT%"
set "EXIT_CODE=%ERRORLEVEL%"

echo.
if "%EXIT_CODE%"=="0" (
    echo Registration batch completed without reported BRAINSFit failures.
) else (
    echo Registration batch finished with exit code %EXIT_CODE%.
    echo Check BRAINSFit_batch_registration_log.csv for details.
)

pause
exit /b %EXIT_CODE%