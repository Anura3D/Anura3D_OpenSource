#      CALCULATE.BAT file
#      print data in the .dat calculation file (instead of a classic .bas template)
proc Anura3D::WriteCalculationFile { filename } {
GiD_WriteCalculationFile init -mode append $filename
set project_path [GiD_Info Project ModelName]
set model_name [file tail $project_path]
set exe_name [GiD_Info Project ProblemType]
GiD_WriteCalculationFile puts -nonewline "\""

# Find the position of the last '/' in the string
set lastSlashIndex [string last / $exe_name]

# Extract the name of the executable after the last '/' symbol
set only_exe_name [string range $exe_name [expr {$lastSlashIndex + 1}] end]

GiD_WriteCalculationFile puts -nonewline $only_exe_name
GiD_WriteCalculationFile puts -nonewline ".exe\" "
GiD_WriteCalculationFile puts -nonewline "\""
GiD_WriteCalculationFile puts -nonewline $project_path
GiD_WriteCalculationFile puts -nonewline ".A3D\\"
GiD_WriteCalculationFile puts -nonewline $model_name
GiD_WriteCalculationFile puts "\""
GiD_WriteCalculationFile puts "PAUSE"
GiD_WriteCalculationFile end
}