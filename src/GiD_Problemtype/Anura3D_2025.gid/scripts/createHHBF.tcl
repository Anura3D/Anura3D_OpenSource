#      Write HHBF file


proc Anura3D::WriteHHBFFile { stageNode projectPath projectName } {            
        
    set xp {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/}
    append xp {container[@n="Hydraulic_head"]/container[@n="Hydraulic_Head_In_Time"]/}
    append xp {container[@n="hydraulic_head_data"]}
       
    set filename [file join [file dir $projectPath] $projectName.A3D $projectName.HHBF]
    if {[file exists $filename]} { file delete -force $filename }           
    set hydraulicheadintimeNode [$stageNode selectNodes $xp]
    if { $hydraulicheadintimeNode == "" } { return }
    GiD_WriteCalculationFile init $filename
    foreach valueNode [$hydraulicheadintimeNode selectNodes {value}] {
        GiD_WriteCalculationFile puts [$valueNode @v]
    }
    GiD_WriteCalculationFile end
         
}
