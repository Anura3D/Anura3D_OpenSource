#      Write PSF files


proc Anura3D::WritePSFFiles { stageNode projectPath projectName } {            
    
    
    set xp {container[@n="Initial_cond"]/container[@n="Stress_initialization"]/}
    append xp {container[@n="Not-horizontal"]/container[@n="Phreatic_surface"]/}
    append xp {blockdata[@n="Phreatic_surface"]/container[@n="water_table_data"]}
    
    # Delete all stored PSF_$icount files
    foreach file [glob -nocomplain -directory [file join [file dir $projectPath] $projectName.A3D] *.PSF_*] {            
        if {[file exists $file]} { file delete -force $file }        
    }
    
    set icount 1
    foreach watertabledataNode [$stageNode selectNodes $xp] {
        set filename [file join [file dir $projectPath] $projectName.A3D $projectName.PSF_$icount]       
        GiD_WriteCalculationFile init $filename
        set numcoord [llength [$watertabledataNode selectNodes {value}]]
        GiD_WriteCalculationFile puts $numcoord
        foreach valueNode [$watertabledataNode selectNodes {value}] {
            GiD_WriteCalculationFile puts [$valueNode @v]
        }
        GiD_WriteCalculationFile end
        incr icount
    }      
}
