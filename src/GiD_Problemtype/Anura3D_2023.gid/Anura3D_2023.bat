@ECHO OFF
REM create input files for A3D calculation
copy %1.dat calculate.bat
REM erase %1.dat
copy %1-1.dat %1.CPS_001
REM erase %1-1.dat
copy %1-2.dat %1.GOM
REM erase %1-2.dat
copy %1-3.dat %1.OPD
REM erase %1-3.dat
REM create A3D folder and copy input files, executable and dll's into it
md ..\%1.A3D
copy calculate.bat "..\%1.A3D\calculate.bat"
copy %1.CPS_001 "..\%1.A3D\%1.CPS_001"
copy %1.GOM "..\%1.A3D\%1.GOM"
REM copy %1.OPD "..\%1.A3D\%1.OPD"
copy %3\exec\Anura3D_2023.exe "..\%1.A3D\Anura3D_2023.exe"
xcopy %3\dll\*.dll "..\%1.A3D\*.dll" /sy