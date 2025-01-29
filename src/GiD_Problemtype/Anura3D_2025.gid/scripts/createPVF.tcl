#      Write PVF file

proc Anura3D::WritePVFFile { stageNode projectPath projectName } {            
        
    set xp {container[@n="BC"]/container[@n="Prescribed_velocity"]/}
    append xp {condition[@n="Nodal_prescribed_velocity"]/group}
    
    foreach groupNode [$stageNode selectNodes $xp] {
        set name [$stageNode @name]
        set n [$groupNode @n]
        
        if {[[$groupNode selectNodes {value[@n="X_velocity_[m/s]"]}] @v] eq "1.0"} {
            set xp {value[@n="X_velocity_[m/s]"]/container[@n="velocity_data"]}
            set velocity_dataNode [$groupNode selectNodes $xp]
                
            set filename [file join [file dir $projectPath] $projectName.A3D $projectName.PVF_${name}_${n}_x]
            if {[file exists $filename]} { file delete -force $filename }                       
            GiD_WriteCalculationFile init $filename
            set numcoord [llength [$velocity_dataNode selectNodes {value}]]
            GiD_WriteCalculationFile puts $numcoord
            foreach valueNode [$velocity_dataNode selectNodes {value}] {
                GiD_WriteCalculationFile puts [$valueNode @v]
            }
            GiD_WriteCalculationFile end
        }
        if {[[$groupNode selectNodes {value[@n="Y_velocity_[m/s]"]}] @v] eq "1.0"} {
            set xp {value[@n="Y_velocity_[m/s]"]/container[@n="velocity_data"]}
            set velocity_dataNode [$groupNode selectNodes $xp]
                
            set filename [file join [file dir $projectPath] $projectName.A3D $projectName.PVF_${name}_${n}_y]
            if {[file exists $filename]} { file delete -force $filename }                       
            GiD_WriteCalculationFile init $filename
            set numcoord [llength [$velocity_dataNode selectNodes {value}]]
            GiD_WriteCalculationFile puts $numcoord
            foreach valueNode [$velocity_dataNode selectNodes {value}] {
                GiD_WriteCalculationFile puts [$valueNode @v]
            }
            GiD_WriteCalculationFile end
            
        }
        if {[[$groupNode selectNodes {value[@n="Z_velocity_[m/s]"]}] @v] eq "1.0"} {
            set xp {value[@n="Z_velocity_[m/s]"]/container[@n="velocity_data"]}
            set velocity_dataNode [$groupNode selectNodes $xp]
                
            set filename [file join [file dir $projectPath] $projectName.A3D $projectName.PVF_${name}_${n}_z]
            if {[file exists $filename]} { file delete -force $filename }                       
            GiD_WriteCalculationFile init $filename
            set numcoord [llength [$velocity_dataNode selectNodes {value}]]
            GiD_WriteCalculationFile puts $numcoord
            foreach valueNode [$velocity_dataNode selectNodes {value}] {
                GiD_WriteCalculationFile puts [$valueNode @v]
            }
            GiD_WriteCalculationFile end
        }
    }              
}
