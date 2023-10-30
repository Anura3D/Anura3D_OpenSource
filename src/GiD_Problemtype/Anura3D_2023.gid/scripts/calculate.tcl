#      CALCULATE.BAT file
#      print data in the .dat calculation file (instead of a classic .bas template)
#      removed the full project path and replaced with the relative path to the current folder of the "calculate.bat" file. 
proc Anura3D::WriteCalculationFile { filename } {
GiD_WriteCalculationFile init -mode append $filename
set project_path [GiD_Info Project ModelName]
set model_name [file tail $project_path]
set exe_name [GiD_Info Project ProblemType]
GiD_WriteCalculationFile puts -nonewline "\".\\"
GiD_WriteCalculationFile puts -nonewline $exe_name
GiD_WriteCalculationFile puts -nonewline ".exe\" "
GiD_WriteCalculationFile puts -nonewline "\".\\"
GiD_WriteCalculationFile puts -nonewline $model_name
GiD_WriteCalculationFile puts "\""
GiD_WriteCalculationFile puts "PAUSE"
GiD_WriteCalculationFile end
}
