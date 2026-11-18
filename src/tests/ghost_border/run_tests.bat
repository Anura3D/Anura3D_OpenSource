@echo off
rem Build and run the ModGhostBorder unit harness with ifx (Debug kernel flags: real kind 8, fpe0, checks).
rem Usage: run_tests.bat <build folder>
setlocal
set OUT=%~1
if "%OUT%"=="" set OUT=%TEMP%\ghost_border_tests
if not exist "%OUT%" mkdir "%OUT%"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" intel64 >nul
set SRC=%~dp0..\..
pushd "%OUT%"
ifx /nologo /free /fpp /real-size:64 /fpe:0 /fp:strict /check:all /traceback "%SRC%\GlobalConstants.FOR" "%SRC%\GhostBorderDetection.FOR" "%~dp0test_ghost_border.f90" /exe:test_ghost_border.exe
if errorlevel 1 (popd & echo BUILD FAILED & exit /b 2)
"%OUT%\test_ghost_border.exe"
set RC=%errorlevel%
popd
exit /b %RC%
