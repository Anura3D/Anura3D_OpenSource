@echo off
LIB /OUT:"Debug\DynamicMPM.lib"  Debug\*.obj
if errorlevel 1 goto IF_BAT_Error
goto IF_BAT_End
:IF_BAT_Error
echo Project : error PRJ0019: A tool returned an error code from "Performing Pre-Link Event..."
exit 1
:IF_BAT_End
