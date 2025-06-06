#      GOM file
#      print data in the .dat calculation file (instead of a classic .bas template)


proc Anura3D::WriteCalculationFile_GOM { filename stageNode project_path model_name } {
    variable current_xml_root
    
    GiD_WriteCalculationFile init $filename
    set exe_name [GiD_Info Project ProblemType]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data
    
    set current_xml_root $root
    
    GiD_WriteCalculationFile puts "### Anura3D_2025 ###"
    
    # DIMENSION
    GiD_WriteCalculationFile puts {$$DIMENSION}
    set dim_path {string(container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set dim_type [$root selectNodes $dim_path]
    if {$dim_type == "2D:plane-strain"} {
        GiD_WriteCalculationFile puts "2D-plane_strain"
    } elseif {$dim_type == "2D:Axissymmetric"} {
        GiD_WriteCalculationFile puts "2D-axisymmetric"
    } elseif {$dim_type == "3D"} {
        GiD_WriteCalculationFile puts "3D-cartesian"
    } elseif {$dim_type == "3D:Axissymmetric"} {
        GiD_WriteCalculationFile puts "3D-cylindrical"
    }
    
    # ELEMENT TYPE
    GiD_WriteCalculationFile puts {$$ELEMENTTYPE}
    set info_mesh [GiD_Mesh get element 1]
    set elem_type [lindex $info_mesh  1]
    set num_nodes [lindex $info_mesh  2]
    if { ($elem_type == "Triangle") && ($num_nodes == "3") } {
        GiD_WriteCalculationFile puts "triangular_3-noded"
    } elseif { ($elem_type == "Tetrahedra") && ($num_nodes == "10") } {
        GiD_WriteCalculationFile puts "tetrahedral_old"
    } else {error [= "INPUT ERROR: Element type not properly defined. Only the following element types are supported: triangular 3-noded, tetrahedral 10-noded."]}   
    
    # FORMULATION
    GiD_WriteCalculationFile puts {$$FORMULATION}
    set layer_path {string(container[@n="General_data"]/value[@n="NLAYERS"]/@v)}
    set layer_type [$root selectNodes $layer_path]
    if {$layer_type == "Single_point"} {
        GiD_WriteCalculationFile puts "single-point"
    } elseif {$layer_type == "Double_point"} {
        GiD_WriteCalculationFile puts "double-point"
    }
    
    # COUNTERS
    set num_nodes [GiD_Info Mesh NumNodes]
    set num_elements [GiD_Info Mesh NumElements]
    GiD_WriteCalculationFile puts {$$STARTCOUNTERS}
    GiD_WriteCalculationFile puts [= "%s %s" $num_elements $num_nodes]

    # NODAL COORDINATES
    set Nodes [GiD_Info Mesh nodes -sublist]
    GiD_WriteCalculationFile puts {$$STARTNODES}
    if {$dim_type == "2D:plane-strain"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $Nodes  $i 1]
            set ycoor [lindex $Nodes  $i 2]
            GiD_WriteCalculationFile puts [= "%s %s" $xcoor $ycoor]}
    } elseif {$dim_type == "2D:Axissymmetric"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $Nodes  $i 1]
            set ycoor [lindex $Nodes  $i 2]
            GiD_WriteCalculationFile puts [= "%s %s" $xcoor $ycoor]}
    } elseif {$dim_type == "3D"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $Nodes  $i 1]
            set ycoor [lindex $Nodes  $i 2]
            set zcoor [lindex $Nodes  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]}
    } elseif {$dim_type == "3D:Axissymmetric"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $Nodes  $i 1]
            set ycoor [lindex $Nodes  $i 2]
            set zcoor [lindex $Nodes  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]}}
    
    # ELEMENT CONNECTIVITIES
    GiD_WriteCalculationFile puts {$$STARTELEMCON}
    if {$dim_type == "2D:plane-strain"} {
        set ElementList [GiD_Info Mesh Elements Triangle -sublist]
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $ElementList  $i 1]
            set ycoor [lindex $ElementList  $i 2]
            set zcoor [lindex $ElementList  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]}
    } elseif {$dim_type == "2D:Axissymmetric"} {
        set ElementList [GiD_Info Mesh Elements Triangle -sublist]
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $ElementList  $i 1]
            set ycoor [lindex $ElementList  $i 2]
            set zcoor [lindex $ElementList  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]}
    } elseif {$dim_type == "3D"} {
        set ElementList [GiD_Info Mesh Elements Tetrahedra -sublist]
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $ElementList  $i 1]
            set ycoor [lindex $ElementList  $i 2]
            set zcoor [lindex $ElementList  $i 3]
            set wcoor [lindex $ElementList  $i 4]
            set jcoor [lindex $ElementList  $i 5]
            set qcoor [lindex $ElementList  $i 6]
            set kcoor [lindex $ElementList  $i 7]
            set lcoor [lindex $ElementList  $i 8]
            set mcoor [lindex $ElementList  $i 9]
            set ncoor [lindex $ElementList  $i 10]
            GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s %s %s %s %s" $xcoor $ycoor $zcoor $wcoor $jcoor $qcoor $kcoor $lcoor $mcoor $ncoor]}
    } elseif {$dim_type == "3D:Axissymmetric"} {
        set ElementList [GiD_Info Mesh Elements Tetrahedra -sublist]
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $ElementList  $i 1]
            set ycoor [lindex $ElementList  $i 2]
            set zcoor [lindex $ElementList  $i 3]
            set wcoor [lindex $ElementList  $i 4]
            set jcoor [lindex $ElementList  $i 5]
            set qcoor [lindex $ElementList  $i 6]
            set kcoor [lindex $ElementList  $i 7]
            set lcoor [lindex $ElementList  $i 8]
            set mcoor [lindex $ElementList  $i 9]
            set ncoor [lindex $ElementList  $i 10]
            GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s %s %s %s %s" $xcoor $ycoor $zcoor $wcoor $jcoor $qcoor $kcoor $lcoor $mcoor $ncoor]}}
    
    # FIXITIES
    # Surface conditions
    ### Solid
    set ov_type "surface"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="solid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_SOLID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    ### Liquid
    set ov_type "surface"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="liquid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_LIQUID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    ### Gas
    set ov_type "surface"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="gas_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_GAS}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats 
    
    # Line conditions
    ### Solid
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="solid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {      
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }        
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_SOLID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    ### Liquid
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="liquid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {      
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }        
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_LIQUID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats 
    ### Gas
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="gas_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {      
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }        
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_GAS}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    
    # Point conditions
    ### Solid
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="solid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_SOLID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    ### Liquid
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="liquid_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_LIQUID}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats       
    ### Gas
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="gas_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_GAS}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile nodes $formats
    
    # NODAL PRESCRIBED VELOCITY (2D/3D)
    # PVF
    # On point
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_velocity_node"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"}
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_POINT}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats
    }
    
    # On line
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_velocity_node"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"}
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_LINE}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats
    }
    
    # On surface
    set ov_type "surface"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_velocity_node"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"}
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_SURFACE}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats 
    }
    
    # On volume 
    set ov_type "volume"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
        set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
        set v3 [$gNode selectNodes {string(value[@n="Z_velocity_node"]/@v)}]
        set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
        set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
        set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
        dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_VOLUME}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats
    }
    
    # MATERIAL POINTS PRESCRIBED VELOCITY 
    # 2D - On surface
    set ov_type "surface"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="MP_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        set v1 [$gNode selectNodes {string(value[@n="X_velocity_MP"]/@v)}]
        set v2 [$gNode selectNodes {string(value[@n="Y_velocity_MP"]/@v)}]
        set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
        set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
        dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile elements -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile elements $formats 
    }
    
    # 3D - On volume 
    set ov_type "volume"
    set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="MP_prescribed_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        set v1 [$gNode selectNodes {string(value[@n="X_velocity_MP"]/@v)}]
        set v2 [$gNode selectNodes {string(value[@n="Y_velocity_MP"]/@v)}]
        set v3 [$gNode selectNodes {string(value[@n="Z_velocity_MP"]/@v)}]
        set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
        set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
        set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
        dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile elements -count $formats]
        GiD_WriteCalculationFile puts {$$PRESCRIBED_MATERIAL_POINT_VELOCITY_VOLUME}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile elements $formats
    }
    
    # INITIAL VELOCITY ON MATERIAL 2D/3D
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {container[@n="Initial_cond"]/condition[@n="Initial_MP_velocity"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {       
            set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_direction"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
    }
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile elements -count $formats]
        GiD_WriteCalculationFile puts {$$INITIAL_VELOCITY_MATERIAL_POINT}    
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile elements $formats 
    }
    
    # INITIAL PHREATIC SURFACE 2D/3D
    # From file        
    set xp {container[@n="Initial_cond"]/container[@n="Stress_initialization"]/}
    append xp {container[@n="Not-horizontal"]/}
    append xp {container[@n="Phreatic_surface"]/blockdata}
    set num_tot 0
    
    foreach gNode [$stageNode selectNodes $xp] { 
        set xp_num_mat {string(value[@n="Number_of_materials"]/@v)}
        set num [$gNode selectNodes $xp_num_mat]
        if {$num == "please specify"} {
            set num 0
        } else {
            set num_tot [expr $num + $num_tot] 
        }
    }
    
    if {$num_tot != 0} {
        GiD_WriteCalculationFile puts {$$INITIAL_WATER_SURFACE_FROM_FILE}
        GiD_WriteCalculationFile puts $num_tot
    }
    
    foreach gNode [$stageNode selectNodes $xp] {
        set name [$gNode @name]
        set water_table_data [$gNode selectNodes {container[@n="water_table_data"]}]
        if {$water_table_data ne ""} {           
            if {$name == "Add water table 1"} {
                set type_number 1
            } elseif {$name == "Add water table 2"} {        
                set type_number 2
            } elseif {$name == "Add water table 3"} {
                set type_number 3
            }      
            set num [$gNode selectNodes {string(value[@n="Number_of_materials"]/@v)}]     
            set type_name [$gNode selectNodes {string(value[@n="material_phreatic_surface_file1"]/@v)}]      
            set MATERIAL_ID [find_material_id $type_name $stageNode]
            GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
            if {$num >= "2"} {
                set type_name [$gNode selectNodes {string(value[@n="material_phreatic_surface_file2"]/@v)}]        
                set MATERIAL_ID [find_material_id $type_name $stageNode]
                GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
            }        
            if {$num >= "3"} {
                set type_name [$gNode selectNodes {string(value[@n="material_phreatic_surface_file3"]/@v)}]
                set MATERIAL_ID [find_material_id $type_name $stageNode]
                GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
            } 
            if {$num >= "4"} {
                set type_name [$gNode selectNodes {string(value[@n="material_phreatic_surface_file4"]/@v)}]
                set MATERIAL_ID [find_material_id $type_name $stageNode]
                GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
            } 
            if {$num == "5"} { 
                set type_name [$gNode selectNodes {string(value[@n="material_phreatic_surface_file5"]/@v)}]
                set MATERIAL_ID [find_material_id $type_name $stageNode]
                GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]       
            }                
        }              
    }
    
    # CONTACT PROPERTIES 2D/3D
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {container[@n="Contact"]/condition[@n="Contact_properties"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {       
        set v1 [$gNode selectNodes {string(value[@n="Number_of_materials"]/@v)}]
        set v2 [$gNode selectNodes {string(value[@n="MATERIAL_1"]/@v)}]
        set v3 [$gNode selectNodes {string(value[@n="Friction_1"]/@v)}]
        set v4 [$gNode selectNodes {string(value[@n="Adhesion_1"]/@v)}]
        set v5 [$gNode selectNodes {string(value[@n="MATERIAL_2"]/@v)}]
        set v6 [$gNode selectNodes {string(value[@n="Friction_2"]/@v)}]
        set v7 [$gNode selectNodes {string(value[@n="Adhesion_2"]/@v)}]
        set v8 [$gNode selectNodes {string(value[@n="MATERIAL_3"]/@v)}]
        set v9 [$gNode selectNodes {string(value[@n="Friction_3"]/@v)}]
        set v10 [$gNode selectNodes {string(value[@n="Adhesion_3"]/@v)}]
        set v11 [$gNode selectNodes {string(value[@n="MATERIAL_4"]/@v)}]
        set v12 [$gNode selectNodes {string(value[@n="Friction_4"]/@v)}]
        set v13 [$gNode selectNodes {string(value[@n="Adhesion_4"]/@v)}]
        if {$v1 == 1} {
            set v5 "NAN"
            set v8 "NAN"
            set v11 "NAN" 
        } elseif {$v1 == 2} {
            set v8 "NAN"
            set v11 "NAN"
        } elseif {$v1 == 3} {
            set v11 "NAN"
        }            
        dict set formats [$gNode @n] "%d $v1 \"$v2\" $v3 $v4 \"$v5\" $v6 $v7 \"$v8\" $v9 $v10 \"$v11\" $v12 $v13\n"}
    
    if { [dict size $formats] } {
        set num [GiD_WriteCalculationFile elements -count $formats]
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            GiD_WriteCalculationFile puts {$$START_BODY_CONTACT_2D}
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} { 
            GiD_WriteCalculationFile puts {$$START_CONTACT_VOLUME}
        }   
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile elements $formats 
    }
    
    # EXCAVATION
    set excavation 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {condition[@n="Solid_Excavation"]/group[@ov=%s]} $ov_type]
    foreach gNode [$stageNode selectNodes $xp] {
        set excavation 1}
    if {$excavation == 1} {
        GiD_WriteCalculationFile puts {$$START_EXCAVATION_SOLID}
    }
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
        set gEntities_id [GiD_EntitiesGroups get $gName $ov_type]
        for {set i 0} {$i < $gEntities_num } {incr i} {
            set id_entity [lindex $gEntities_id $i]
            set elements [GiD_Geometry get $ov_type $id_entity mesh]
            set elements_id [lindex $elements 4]
            set elements_num [llength $elements_id]
            for {set j 0} {$j < $elements_num } {incr j} {
                set id_element [lindex $elements_id $j]
                GiD_WriteCalculationFile puts [= "%s %s" $id_entity $id_element]
            } 
        }
    }
    
    # CONSTRUCTION
    set construction 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {container[@n="Construction"]/condition[@n="Solid_Construction"]/group[@ov=%s]} $ov_type]
    foreach gNode [$stageNode selectNodes $xp] {
        set construction 1}
    if {$construction == 1} {
        GiD_WriteCalculationFile puts {$$START_CONSTRUCTION_SOLID}
    }
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
        set gEntities_id [GiD_EntitiesGroups get $gName $ov_type]
        for {set i 0} {$i < $gEntities_num } {incr i} {
            set id_entity [lindex $gEntities_id $i]
            set elements [GiD_Geometry get $ov_type $id_entity mesh]
            set elements_id [lindex $elements 4]
            set elements_num [llength $elements_id]
            for {set j 0} {$j < $elements_num } {incr j} {
                set id_element [lindex $elements_id $j]
                GiD_WriteCalculationFile puts [= "%s %s" $id_entity $id_element]} }
    }
    
    # HYDRAULIC BOUNDARY CONDITIONS
    # Hydraulic head  
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/}  
        append xp_min {value[@n="Minimum_coordinates_hydraulic_head_2D"]}         
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin  
        set zmin 0.0
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/}               
        append xp_max {value[@n="Maximum_coordinates_hydraulic_head_2D"]}
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax  
        set zmin 0.0
        set xp_hydraulic_head {container[@n="Calculation_Data"]/value[@n="HYDRAULIC_HEAD"]}
        set vNode [$stageNode selectNodes $xp_hydraulic_head]
        set v [get_domnode_attribute $vNode v]       
        if { $v } {
            GiD_WriteCalculationFile puts {$$BOUNDARY_HYDRAULIC_HEAD_AREA}                
            GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
            GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
        }
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/}
        append xp_min {value[@n="Minimum_coordinates_hydraulic_head_3D"]}        
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin zmin          
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/}
        append xp_max {value[@n="Maximum_coordinates_hydraulic_head_3D"]}       
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax zmax  
        set xp_hydraulic_head {container[@n="Calculation_Data"]/value[@n="HYDRAULIC_HEAD"]}
        set vNode [$stageNode selectNodes $xp_hydraulic_head]
        set v [get_domnode_attribute $vNode v]
        if { $v } {                
            GiD_WriteCalculationFile puts {$$BOUNDARY_HYDRAULIC_HEAD_AREA}          
            GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
            GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
            GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax] 
        }
    }    
    
    # Seepage face
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/}  
        append xp_min {value[@n="Minimum_coordinates_seepage_face_2D"]}         
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin  
        set zmin 0.0
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/}               
        append xp_max {value[@n="Maximum_coordinates_seepage_face_2D"]}
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax  
        set zmin 0.0
        set xp_seepage_face {container[@n="Calculation_Data"]/value[@n="SEEPAGE_FACE"]}
        set vNode [$stageNode selectNodes $xp_seepage_face]
        set v [get_domnode_attribute $vNode v]
        if { $v } {       
            GiD_WriteCalculationFile puts {$$BOUNDARY_SEEPAGE_AREA}
            GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
            GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
        }       
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/}
        append xp_min {value[@n="Minimum_coordinates_seepage_face_3D"]}        
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin zmin          
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/}
        append xp_max {value[@n="Maximum_coordinates_seepage_face_3D"]}       
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax zmax  
        set xp_seepage_face {container[@n="Calculation_Data"]/value[@n="SEEPAGE_FACE"]}
        set vNode [$stageNode selectNodes $xp_seepage_face]
        set v [get_domnode_attribute $vNode v]        
        if { $v } {
            GiD_WriteCalculationFile puts {$$BOUNDARY_SEEPAGE_AREA}  
            GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
            GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
            GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax]
        }       
    }    
    
    # Infiltration
    set xp [format_xpath {string(container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="X_direction"]/@v)}]
    set x_value [$stageNode selectNodes $xp]
    set yp [format_xpath {string(container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="Y_direction"]/@v)}]
    set y_value [$stageNode selectNodes $yp]
    set zp [format_xpath {string(container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="Z_direction"]/@v)}]  
    set z_value [$stageNode selectNodes $zp]   
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/}  
        append xp_min {value[@n="Minimum_coordinates_infiltration_2D"]}         
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin  
        set zmin 0.0
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/}               
        append xp_max {value[@n="Maximum_coordinates_infiltration_2D"]}
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax  
        set zmin 0.0
        set xp_infiltration {container[@n="Calculation_Data"]/value[@n="INFILTRATION"]}        
        set vNode [$stageNode selectNodes $xp_infiltration]
        set v [get_domnode_attribute $vNode v]   
        if { $v } {	
        GiD_WriteCalculationFile puts {$$BOUNDARY_INFILTRATION_AREA}		
        GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
        GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
        GiD_WriteCalculationFile puts {$$INFILTRATION_RATE}
        GiD_WriteCalculationFile puts [= "%s %s" $x_value $y_value]
		}
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set xp_min {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/}
        append xp_min {value[@n="Minimum_coordinates_infiltration_3D"]}        
        set vNode [$stageNode selectNodes $xp_min]
        set v [split [$vNode @v] ","]
        lassign $v xmin ymin zmin          
        set xp_max {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/}
        append xp_max {value[@n="Maximum_coordinates_infiltration_3D"]}       
        set vNode [$stageNode selectNodes $xp_max]
        set v [split [$vNode @v] ","]
        lassign $v xmax ymax zmax  
        set xp_infiltration {container[@n="Calculation_Data"]/value[@n="INFILTRATION"]}
        set vNode [$stageNode selectNodes $xp_infiltration]
        set v [get_domnode_attribute $vNode v]
        if { $v } {		
        GiD_WriteCalculationFile puts {$$BOUNDARY_INFILTRATION_AREA}		
        GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
        GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
        GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax] 
        GiD_WriteCalculationFile puts {$$INFILTRATION_RATE}
        GiD_WriteCalculationFile puts [= "%s %s" $x_value $y_value $z_value]
		}
    }  
    
    # Reaction forces 2D/3D  (2D/3D)
    # 2D On line
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "line"
        set xp [format_xpath {condition[@n="Reaction_forces"]/group[@ov=%s]} $ov_type]
        set formats ""
        foreach gNode [$stageNode selectNodes $xp] {
            set line_name [$gNode selectNodes {string(value[@n="line_identifier"]/@v)}]
            dict set formats [$gNode @n] "\"$line_name\" %d %d %d\n"
        }    
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "surface"
        set xp [format_xpath {condition[@n="Reaction_forces"]/group[@ov=%s]} $ov_type]
        set formats ""
        foreach gNode [$stageNode selectNodes $xp] {
            set line_name [$gNode selectNodes {string(value[@n="line_identifier"]/@v)}]
            dict set formats [$gNode @n] "\"$line_name\" %d %d %d %d %d %d %d\n"
        }
    }
    
    set num [GiD_WriteCalculationFile elements -count $formats]
    GiD_WriteCalculationFile puts {$$START_OUTPUT_REACTION_FORCES}    
    GiD_WriteCalculationFile puts $num
    GiD_WriteCalculationFile elements -print_faces_conecs $formats
    
    
    # ABSORBING BOUNDARIES 2D/3D
    # Surface conditions
    ### Solid
    set num1 0
    if {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="solid_absorbing_surface"]/group[@ov=%s]} $ov_type]
        set formats ""
        foreach gNode [$stageNode selectNodes $xp] {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
        }
        set num1 [GiD_WriteCalculationFile nodes -count $formats]
        if {$num1 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_SURFACE_SOLID}    
            GiD_WriteCalculationFile puts $num1
            GiD_WriteCalculationFile nodes $formats
        }
    }
    # line conditions
    ### Solid
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="solid_absorbing_line"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"}
    }
    set num2 [GiD_WriteCalculationFile nodes -count $formats]
    if {$num2 != 0} {
        GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_LINE_SOLID}    
        GiD_WriteCalculationFile puts $num2
        GiD_WriteCalculationFile nodes $formats
    }
    # node conditions
    ### Solid
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="solid_absorbing_node"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"}
    }
    set num3 [GiD_WriteCalculationFile nodes -count $formats]
    if {$num3 != 0} {
        GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_POINT_SOLID}    
        GiD_WriteCalculationFile puts $num3
        GiD_WriteCalculationFile nodes $formats
    }
    
    # Surface conditions
    ### Liquid
    set num4 0
    if {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {        
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="fluid_absorbing_surface"]/group[@ov=%s]} $ov_type]
        set formats ""
        foreach gNode [$stageNode selectNodes $xp] {
            
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_l"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha_l"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta_l"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint_l"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha_l"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta_l"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint_l"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha_l"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta_l"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
        }
        set num4 [GiD_WriteCalculationFile nodes -count $formats]
        if {$num4 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_SURFACE_LIQUID}    
            GiD_WriteCalculationFile puts $num4
            GiD_WriteCalculationFile nodes $formats
        }
    }
    # line conditions
    ### Solid
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="fluid_absorbing_line"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_l"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha_l"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta_l"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint_l"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha_l"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta_l"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_l"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha_l"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta_l"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint_l"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha_l"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta_l"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint_l"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha_l"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta_l"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"}
    }
    
    set num5 [GiD_WriteCalculationFile nodes -count $formats]
    if {$num5 != 0} {
        GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_LINE_LIQUID}    
        GiD_WriteCalculationFile puts $num5
        GiD_WriteCalculationFile nodes $formats
    }
    
    # node conditions
    ### Solid
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/condition[@n="fluid_absorbing_node"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_l"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha_l"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta_l"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint_l"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha_l"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta_l"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
            set v1 [$gNode selectNodes {string(value[@n="X_Constraint_l"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="X_alpha_l"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_delta_l"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_Constraint_l"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_alpha_l"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Y_delta_l"]/@v)}]
            set v7 [$gNode selectNodes {string(value[@n="Z_Constraint_l"]/@v)}]
            set v8 [$gNode selectNodes {string(value[@n="Z_alpha_l"]/@v)}]
            set v9 [$gNode selectNodes {string(value[@n="Z_delta_l"]/@v)}]
            dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"}
    }
    set num6 [GiD_WriteCalculationFile nodes -count $formats]
    if {$num6 != 0} {
        GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_POINT_LIQUID}    
        GiD_WriteCalculationFile puts $num6
        GiD_WriteCalculationFile nodes $formats
    }
    if {$num1 !=0 || $num2 !=0 || $num3 !=0 || $num4 !=0 || $num5 !=0 || $num6 !=0} {
        set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/value[@n="material_absorbing"]}]
        set node [$stageNode selectNodes $xp]
        set material_name [$node getAttribute "v"]
        set MATERIAL_ID [find_material_id $material_name $stageNode]
        GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_REFERENCE_MATERIAL_INDEX}    
        GiD_WriteCalculationFile puts $MATERIAL_ID
    }    
    
    # Loading_Conditions (2D/3D)
    
    # 2D On line
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "line"
        set xp [format_xpath {container[@n="BC"]/container[@n="Loading_Conditions"]/condition[@n="Solid_traction"]/group[@ov=%s]} $ov_type]
        set formats ""
        set formats_b ""
        set formats_mp ""
        set formats_b_mp ""
        foreach gNode [$stageNode selectNodes $xp] {        
            
            set l_dis [$gNode selectNodes {string(value[@n="loading_distribution"]/@v)}]
            set local [$gNode selectNodes {string(value[@n="Local_axes"]/@v)}]
            set load_on [$gNode selectNodes {string(value[@n="apply_traction_on"]/@v)}]
            set load_class [$gNode selectNodes {string(value[@n="load_system_"]/@v)}]
            
            set list_group [$gNode @n]
            set nodes [GiD_EntitiesGroups get $list_group faces]
            set elements_array [lindex $nodes 0]
            set face_array [lindex $nodes 1]
            set num_el [llength $elements_array]
            
            if {$load_on == "nodes"} {
                if {$l_dis == "uniform"} {
                    if {$local != "0"} {
                    }
                    set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
                    set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
                    set v3 $v1
                    set v4 $v2
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord"]/@v)}]
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX_of_tractionX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY_of_tractionX"]/@v)}]
                    set gradyy [$gNode selectNodes {string(value[@n="gradientY_of_tractionY"]/@v)}]
                    set gradyx [$gNode selectNodes {string(value[@n="gradientX_of_tractionY"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set v2 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)]
                        set cntwo [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v3 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set v4 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v4 $v1 $v2"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v4 $v1 $v2"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                        
                    } 
                }
                
                
                
            } elseif {$load_on == "material points"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
                    set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
                    set v3 $v1
                    set v4 $v2
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord"]/@v)}]
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX_of_tractionX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY_of_tractionX"]/@v)}]
                    set gradyy [$gNode selectNodes {string(value[@n="gradientY_of_tractionY"]/@v)}]
                    set gradyx [$gNode selectNodes {string(value[@n="gradientX_of_tractionY"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set v2 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v3 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set v4 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v4 $v1 $v2"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v4 $v1 $v2"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                        
                    } 
                }                
                
                
            }
        }
    }        elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {        
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Loading_Conditions"]/condition[@n="Solid_traction"]/group[@ov=%s]} $ov_type]
        set formats ""
        set formats_b ""
        set formats_mp ""
        set formats_b_mp ""
        foreach gNode [$stageNode selectNodes $xp] {       
            
            set l_dis [$gNode selectNodes {string(value[@n="loading_distribution"]/@v)}]
            set local [$gNode selectNodes {string(value[@n="Local_axes"]/@v)}]
            set load_on [$gNode selectNodes {string(value[@n="apply_traction_on"]/@v)}]
            set load_class [$gNode selectNodes {string(value[@n="load_system_"]/@v)}]
            
            set list_group [$gNode @n]
            set nodes [GiD_EntitiesGroups get $list_group faces]
            set elements_array [lindex $nodes 0]
            set face_array [lindex $nodes 1]
            set num_el [llength $elements_array]
            
            if {$load_on == "nodes"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
                    set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
                    set v3 [$gNode selectNodes {string(value[@n="Z_direction"]/@v)}]
                    if {$local != "0"} {
                    }                                                
                    set v4 $v1
                    set v5 $v2
                    set v6 $v3
                    set v7 $v1
                    set v8 $v2
                    set v9 $v3
                    set v10 $v1
                    set v11 $v2
                    set v12 $v3
                    set v13 $v1
                    set v14 $v2
                    set v15 $v3
                    set v16 $v1
                    set v17 $v2
                    set v18 $v3
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord"]/@v)}]
                    set z0 [$gNode selectNodes {string(value[@n="reference_point_Z_coord"]/@v)}]
                    
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set tracZ0 [$gNode selectNodes {string(value[@n="tractionZ_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX_of_tractionX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY_of_tractionX"]/@v)}]
                    set gradyy [$gNode selectNodes {string(value[@n="gradientY_of_tractionY"]/@v)}]
                    set gradyx [$gNode selectNodes {string(value[@n="gradientX_of_tractionY"]/@v)}]
                    set gradxz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionX"]/@v)}]
                    set gradyz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionY"]/@v)}]
                    set gradzx [$gNode selectNodes {string(value[@n="gradientX_of_tractionZ"]/@v)}]
                    set gradzy [$gNode selectNodes {string(value[@n="gradientY_of_tractionZ"]/@v)}]
                    set gradzz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionZ"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v2 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v3 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v4 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v5 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v6 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 2] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v7 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v8 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v9 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 3] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v10 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v11 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v12 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 4] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v13 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v14 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v15 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 5] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v16 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v17 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v18 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]                                                                                                
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }                                                                                                
                        
                    } 
                }
                
                
                
            } elseif {$load_on == "material points"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="X_direction"]/@v)}]
                    set v2 [$gNode selectNodes {string(value[@n="Y_direction"]/@v)}]
                    set v3 [$gNode selectNodes {string(value[@n="Z_direction"]/@v)}]
                    set v4 $v1
                    set v5 $v2
                    set v6 $v3
                    set v7 $v1
                    set v8 $v2
                    set v9 $v3
                    set v10 $v1
                    set v11 $v2
                    set v12 $v3
                    set v13 $v1
                    set v14 $v2
                    set v15 $v3
                    set v16 $v1
                    set v17 $v2
                    set v18 $v3
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord"]/@v)}]
                    set z0 [$gNode selectNodes {string(value[@n="reference_point_Z_coord"]/@v)}]
                    
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set tracZ0 [$gNode selectNodes {string(value[@n="tractionZ_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX_of_tractionX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY_of_tractionX"]/@v)}]
                    set gradyy [$gNode selectNodes {string(value[@n="gradientY_of_tractionY"]/@v)}]
                    set gradyx [$gNode selectNodes {string(value[@n="gradientX_of_tractionY"]/@v)}]
                    set gradxz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionX"]/@v)}]
                    set gradyz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionY"]/@v)}]
                    set gradzx [$gNode selectNodes {string(value[@n="gradientX_of_tractionZ"]/@v)}]
                    set gradzy [$gNode selectNodes {string(value[@n="gradientY_of_tractionZ"]/@v)}]
                    set gradzz [$gNode selectNodes {string(value[@n="gradientZ_of_tractionZ"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v2 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v3 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v4 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v5 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v6 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 2] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v7 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v8 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v9 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 3] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v10 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v11 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v12 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 4] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v13 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v14 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v15 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 5] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v16 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v17 [expr $tracY0+$gradyx*($X1-$x0)+$gradyy*($Y1-$y0)+$gradyz*($Z1-$z0)]
                        set v18 [expr $tracY0+$gradzx*($X1-$x0)+$gradzy*($Y1-$y0)+$gradzz*($Z1-$z0)]                                                                                                
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }                                                                                                
                        
                    } 
                }
                
                
            }
        }
    }
    
    set num [llength $formats]
    if { $num != 0 } {
        set num [llength $formats]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_NODES_SOLID}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats $i]
        }
    }
    set num [llength $formats_b]
    if { $num != 0 } {
        set num [llength $formats_b]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_NODES_SOLID_B}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_b $i]
        }
    }
    
    set num [llength $formats_mp]
    if { $num != 0 } {
        set num [llength $formats_mp]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_MATERIAL_POINTS_SOLID}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_mp $i]
        }
    }
    
    set num [llength $formats_b_mp]
    if { $num != 0 } {
        set num [llength $formats_b_mp]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_MATERIAL_POINTS_SOLID_B}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_b_mp $i]
        }
    }
    
    
    
    
    # Loading_Conditions (2D/3D)
    # Pressure
    # 2D On line
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "line"
        set xp [format_xpath {container[@n="BC"]/container[@n="Loading_Conditions"]/condition[@n="Liquid_Pressure"]/group[@ov=%s]} $ov_type]
        set formats ""
        set formats_b ""
        set formats_mp ""
        set formats_b_mp ""
        foreach gNode [$stageNode selectNodes $xp] {        
            
            set l_dis [$gNode selectNodes {string(value[@n="all_directions"]/@v)}]
            set load_on [$gNode selectNodes {string(value[@n="apply_pressure_on"]/@v)}]
            set load_class [$gNode selectNodes {string(value[@n="load_system"]/@v)}]
            
            set list_group [$gNode @n]
            set nodes [GiD_EntitiesGroups get $list_group faces]
            set elements_array [lindex $nodes 0]
            set face_array [lindex $nodes 1]
            set num_el [llength $elements_array]
            
            if {$load_on == "nodes"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="pressure"]/@v)}]
                    set v2 $v1
                    set v3 $v1
                    set v4 $v1
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord_"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord_"]/@v)}]
                    set tracX0 [$gNode selectNodes {string(value[@n="pressure_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set cntwo [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cntwo 0]
                        set Y1 [lindex $cntwo 1]
                        set v3 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v1 $v3 $v3"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v3 $v1 $v1"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v1 $v3 $v3"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v3 $v1 $v1"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                        
                    } 
                }
                
                
                
            } elseif {$load_on == "material points"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="pressure"]/@v)}]
                    set v2 $v1
                    set v3 $v1
                    set v4 $v1
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v2 $v3 $v4"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v1 $v2 $v3 $v4"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord_"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord_"]/@v)}]
                    set tracX0 [$gNode selectNodes {string(value[@n="pressure_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY"]/@v)}]
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        set cntwo [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cntwo 0]
                        set Y1 [lindex $cntwo 1]
                        set v3 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v1 $v3 $v3"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v3 $v1 $v1"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] $v1 $v1 $v3 $v3"
                            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0] $v3 $v3 $v1 $v1"
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                if {$string_to_test eq [lindex $formats_b_mp $j]} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                        
                    } 
                }              
                
                
            }
        }
    }        elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {        
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Loading_Conditions"]/condition[@n="Liquid_Pressure"]/group[@ov=%s]} $ov_type]
        set formats ""
        set formats_b ""
        set formats_mp ""
        set formats_b_mp ""
        foreach gNode [$stageNode selectNodes $xp] {       
            
            set l_dis [$gNode selectNodes {string(value[@n="all_directions"]/@v)}]
            set load_on [$gNode selectNodes {string(value[@n="apply_pressure_on"]/@v)}]
            set load_class [$gNode selectNodes {string(value[@n="load_system"]/@v)}]
            
            set list_group [$gNode @n]
            set nodes [GiD_EntitiesGroups get $list_group faces]
            set elements_array [lindex $nodes 0]
            set face_array [lindex $nodes 1]
            set num_el [llength $elements_array]
            
            if {$load_on == "nodes"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="pressure"]/@v)}]
                    set v2 $v1
                    set v3 $v1                                             
                    set v4 $v1
                    set v5 $v2
                    set v6 $v3
                    set v7 $v1
                    set v8 $v2
                    set v9 $v3
                    set v10 $v1
                    set v11 $v2
                    set v12 $v3
                    set v13 $v1
                    set v14 $v2
                    set v15 $v3
                    set v16 $v1
                    set v17 $v2
                    set v18 $v3
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord_"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord_"]/@v)}]
                    set z0 [$gNode selectNodes {string(value[@n="reference_point_Z_coord_"]/@v)}]
                    
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set tracZ0 [$gNode selectNodes {string(value[@n="tractionZ_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY"]/@v)}]
                    set gradxz [$gNode selectNodes {string(value[@n="gradientZ"]/@v)}]
                    
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v2 $v1
                        set v3 $v1
                        set cnone [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v4 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v5 $v4
                        set v6 $v4
                        set cnone [GiD_Mesh get node [lindex $nodes_local 2] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v7 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v8 $v7
                        set v9 $v7
                        set cnone [GiD_Mesh get node [lindex $nodes_local 3] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v10 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v11 $v10
                        set v12 $v10
                        set cnone [GiD_Mesh get node [lindex $nodes_local 4] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v13 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v14 $v13
                        set v15 $v13
                        set cnone [GiD_Mesh get node [lindex $nodes_local 5] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v16 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v17 $v16
                        set v18 $v16                                                                                                
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b $currentstring
                            }
                        }                                                                                                
                        
                    } 
                }
                
                
                
            } elseif {$load_on == "material points"} {
                if {$l_dis == "uniform"} {
                    set v1 [$gNode selectNodes {string(value[@n="pressure"]/@v)}]
                    set v2 $v1
                    set v3 $v1                                             
                    set v4 $v1
                    set v5 $v2
                    set v6 $v3
                    set v7 $v1
                    set v8 $v2
                    set v9 $v3
                    set v10 $v1
                    set v11 $v2
                    set v12 $v3
                    set v13 $v1
                    set v14 $v2
                    set v15 $v3
                    set v16 $v1
                    set v17 $v2
                    set v18 $v3
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }
                    }
                } else {
                    
                    
                    set x0 [$gNode selectNodes {string(value[@n="reference_point_X_coord_"]/@v)}]
                    set y0 [$gNode selectNodes {string(value[@n="reference_point_Y_coord_"]/@v)}]
                    set z0 [$gNode selectNodes {string(value[@n="reference_point_Z_coord_"]/@v)}]
                    
                    set tracX0 [$gNode selectNodes {string(value[@n="tractionX_at_reference_point"]/@v)}]
                    set tracY0 [$gNode selectNodes {string(value[@n="tractionY_at_reference_point"]/@v)}]
                    set tracZ0 [$gNode selectNodes {string(value[@n="tractionZ_at_reference_point"]/@v)}]
                    set gradxx [$gNode selectNodes {string(value[@n="gradientX"]/@v)}]
                    set gradxy [$gNode selectNodes {string(value[@n="gradientY"]/@v)}]
                    set gradxz [$gNode selectNodes {string(value[@n="gradientZ"]/@v)}]
                    
                    
                    for {set i 0} {$i < $num_el} {incr i} {
                        set id_el [lindex $elements_array $i]
                        set id_face [lindex $face_array $i]
                        set nodes_local [GiD_Mesh get element $id_el face $id_face]
                        set cnone [GiD_Mesh get node [lindex $nodes_local 0] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v1 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v2 $v1
                        set v3 $v1
                        set cnone [GiD_Mesh get node [lindex $nodes_local 1] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v4 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v5 $v4
                        set v6 $v4
                        set cnone [GiD_Mesh get node [lindex $nodes_local 2] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v7 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v8 $v7
                        set v9 $v7
                        set cnone [GiD_Mesh get node [lindex $nodes_local 3] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v10 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v11 $v10
                        set v12 $v10
                        set cnone [GiD_Mesh get node [lindex $nodes_local 4] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v13 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v14 $v13
                        set v15 $v13
                        set cnone [GiD_Mesh get node [lindex $nodes_local 5] coordinates]
                        set X1 [lindex $cnone 0]
                        set Y1 [lindex $cnone 1]
                        set Z1 [lindex $cnone 2]
                        set v16 [expr $tracX0+$gradxx*($X1-$x0)+$gradxy*($Y1-$y0)+$gradxz*($Z1-$z0)]
                        set v17 $v16
                        set v18 $v16                                                                                                
                        if {$load_class == "A"} {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_mp $currentstring
                            }
                        } else {
                            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1] [lindex $nodes_local 2] [lindex $nodes_local 3] [lindex $nodes_local 4] [lindex $nodes_local 5] $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14 $v15 $v16 $v17 $v18"
                            set nums1 [lsort [lrange [split $currentstring " "] 0 2]]
                            set formatted_nums1 [join $nums1 " "]
                            set flag 0
                            for {set j 0} {$j < $i} {incr j} {
                                set str [lindex $formats_b_mp $j]
                                set nums [lsort [lrange [split $str " "] 0 2]]
                                set formatted_nums2 [join $nums " "]
                                if {$formatted_nums2 eq $formatted_nums1} {
                                    set flag 1
                                }
                            }
                            if {$flag == 0} {
                                lappend formats_b_mp $currentstring
                            }
                        }                                                                                                
                        
                    } 
                }
                
                
            }
        }
    }
    
    set num [llength $formats]
    if { $num != 0 } {
        set num [llength $formats]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_NODES_LIQUID}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats $i]
        }
    }
    set num [llength $formats_b]
    if { $num != 0 } {
        set num [llength $formats_b]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_NODES_LIQUID_B}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_b $i]
        }
    }
    
    set num [llength $formats_mp]
    if { $num != 0 } {
        set num [llength $formats_mp]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_MATERIAL_POINTS_LIQUID}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_mp $i]
        }
    }
    
    set num [llength $formats_b_mp]
    if { $num != 0 } {
        set num [llength $formats_b_mp]
        GiD_WriteCalculationFile puts {$$START_LOAD_ON_MATERIAL_POINTS_LIQUID_B}    
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats_b_mp $i]
        }
    }
    
    # Soil surface (2D)
    # Line conditions
    ### Solid
    set ov_type "line"
    set xp [format_xpath {container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Not-horizontal"]/condition[@n="Soil_surface"]/group[@ov=%s]} $ov_type]
    set formats ""       
    foreach gNode [$stageNode selectNodes $xp] {
        set list_group [$gNode @n]
        set nodes [GiD_EntitiesGroups get $list_group faces]
        set elements_array [lindex $nodes 0]
        set face_array [lindex $nodes 1]
        set num_el [llength $elements_array]
        for {set i 0} {$i < $num_el} {incr i} {
            set id_el [lindex $elements_array $i]
            set id_face [lindex $face_array $i]
            set nodes_local [GiD_Mesh get element $id_el face $id_face]
            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1]"
            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0]"
            set flag 0  
            for {set j 0} {$j < $i} {incr j} {
                if {$string_to_test eq [lindex $formats $j]} {
                    set flag 1}
            }
            if {$flag == 0} {
                lappend formats $currentstring
            }
        }
        set num [llength $formats]
        GiD_WriteCalculationFile puts {$$START_SOIL_SURFACE_NODES}
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats $i]
        }
    }
    
    # Phreatic surface (2D)
    # Line conditions
    ### Solid
    set ov_type "line"
    set xp [format_xpath {container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Not-horizontal"]/container[@n="Phreatic_surface"]/condition[@n="Phreatic_surface_line"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        set list_group [$gNode @n]
        set nodes [GiD_EntitiesGroups get $list_group faces]
        set elements_array [lindex $nodes 0]
        set face_array [lindex $nodes 1]
        set num_el [llength $elements_array]        
        for {set i 0} {$i < $num_el} {incr i} {
            set id_el [lindex $elements_array $i]
            set id_face [lindex $face_array $i]
            set nodes_local [GiD_Mesh get element $id_el face $id_face]
            set currentstring "[lindex $nodes_local 0] [lindex $nodes_local 1]"
            set string_to_test "[lindex $nodes_local 1] [lindex $nodes_local 0]"
            set flag 0  
            for {set j 0} {$j < $i} {incr j} {
                if {$string_to_test eq [lindex $formats $j]} {
                    set flag 1}
            }
            if {$flag == 0} {
                lappend formats $currentstring
            }
        }
        set num [llength $formats]
        GiD_WriteCalculationFile puts {$$START_PHREATIC_SURFACE_NODES}
        GiD_WriteCalculationFile puts $num
        for {set i 0} {$i < $num} {incr i} {
            GiD_WriteCalculationFile puts [lindex $formats $i]
        }
    }
    
    # Moving_mesh
    # Extending_mesh
    set ov_type "point"
    set xp [format_xpath {container[@n="Moving_mesh"]/condition[@n="Extending_mesh"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        dict set formats [$gNode @n] "%d\n"
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    if {$num != 0} {
        GiD_WriteCalculationFile puts {EXTENDING_MESH_CORNER_NODES}    
        GiD_WriteCalculationFile nodes $formats
    }
    # Compressing_mesh
    set ov_type "point"
    set xp [format_xpath {container[@n="Moving_mesh"]/condition[@n="Compressing_mesh"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        dict set formats [$gNode @n] "%d\n"
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    if {$num != 0} {
        GiD_WriteCalculationFile puts {$$COMPRESSING_MESH_CORNER_NODES}    
        GiD_WriteCalculationFile nodes $formats
    }
    # Moving_mesh
    set ov_type "point"
    set xp [format_xpath {container[@n="Moving_mesh"]/condition[@n="Moving_mesh"]/group[@ov=%s]} $ov_type]
    set formats ""
    foreach gNode [$stageNode selectNodes $xp] {
        dict set formats [$gNode @n] "%d\n"
        set node [$gNode selectNodes [format_xpath {value[@n="Mm_dir"]}]]
        set mesh_dir [$node getAttribute "v"]
    }
    set num [GiD_WriteCalculationFile nodes -count $formats]
    if {$num != 0} {
        GiD_WriteCalculationFile puts {$$MOVING_MESH_CORNER_NODES}    
        GiD_WriteCalculationFile nodes $formats
        GiD_WriteCalculationFile puts {$$MOVING_MESH_DIRECTION}    
        GiD_WriteCalculationFile puts $mesh_dir
        
        # Moving_mesh Reference material
        set xp [format_xpath {container[@n="Moving_mesh"]/value[@n="Reference_material"]}]
        set node [$stageNode selectNodes $xp]
        set material_name [$node getAttribute "v"]
        set MATERIAL_ID [find_material_id $material_name $stageNode]
        GiD_WriteCalculationFile puts {$$MOVING_MESH_REFERENCE_MATERIAL_INDEX}    
        GiD_WriteCalculationFile puts $MATERIAL_ID
    }
    
    # MATERIALS
    set xp [format_xpath {container[@n="materials"]/blockdata}]   
    set list_len [llength [$stageNode selectNodes $xp]]
    GiD_WriteCalculationFile puts {$$NUMBER_OF_MATERIALS}    
    GiD_WriteCalculationFile puts $list_len 
    set int 1
    foreach gNode [$stageNode selectNodes $xp] {
        set type [$gNode @name]
        set MATERIAL_ID [find_material_id $type $stageNode]
        GiD_WriteCalculationFile puts {$$MATERIAL_INDEX}
        GiD_WriteCalculationFile puts $MATERIAL_ID                                                 
        GiD_WriteCalculationFile puts {$$MATERIAL_NAME}    
        GiD_WriteCalculationFile puts $type 
        set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="material_type_"]}]]
        set typename [$node getAttribute "v"]
        if {$typename == "Dry material"} {
            set type "dry_material"
        } elseif {$typename == "Saturated material-drained"} {
            set type "saturated_material_drained"
        } elseif {$typename == "Saturated material-undrained effective stress"} {
            set type "saturated_material_undrained_effective"
        } elseif {$typename == "Saturated material-undrained total stress"} {
            set type "saturated_material_undrained_total"
        } elseif {$typename == "Saturated material-fully coupled"} {
            set type "saturated_material_coupled"
        } elseif {$typename == "Unsaturated material-2-phase with suction effect"} {
            set type "unsaturated_material_2phase_suction"
        } elseif {$typename == "Unsaturated material-3-phase fully coupled"} {
            set type "unsaturated_material_2phase_suction"
        } elseif {$typename == "Liquid"} {
            set type "liquid"}
        GiD_WriteCalculationFile puts {$$MATERIAL_TYPE}    
        GiD_WriteCalculationFile puts $type
        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="initial_porosity_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$POROSITY_SOLID}    
            GiD_WriteCalculationFile puts $type
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_solid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_SOLID }    
            GiD_WriteCalculationFile puts $type                                             
        }
        if {$typename != "Dry material"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_LIQUID }    
            GiD_WriteCalculationFile puts $type                                              
        }
        if {$typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_GAS}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {$typename == "Saturated material-fully coupled"||$typename == "Unsaturated material-2-phase with suction effect"||$typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="intrinsic_permeability_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY_LIQUID}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="intrinsic_permeability_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY_GAS}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {($typename == "Liquid") || ($typename == "Saturated material-fully coupled") || ($typename == "Unsaturated material-2-phase with suction effect") || ($typename == "Unsaturated material-3-phase fully coupled")} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="bulk__modulus_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$BULK_MODULUS_LIQUID}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="bulk__modulus_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$BULK_MODULUS_GAS}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {($typename == "Liquid") || ($typename == "Saturated material-fully coupled") || ($typename == "Unsaturated material-2-phase with suction effect") || ($typename == "Unsaturated material-3-phase fully coupled") } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="dinamic_viscosity_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DYNAMIC_VISCOSITY_LIQUID}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="dinamic_viscosity_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DYNAMIC_VISCOSITY_GAS}    
            GiD_WriteCalculationFile puts $type       
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="elastic_swelling_index_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$SWELLING_INDEX}    
            GiD_WriteCalculationFile puts $type    
        }         
        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="K0-value_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$K0_VALUE_SOLID}    
            GiD_WriteCalculationFile puts $type                                                 
        }
        if {$typename == "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="liquid_cavitation_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$LIQUID_CAVITATION}    
            GiD_WriteCalculationFile puts $type
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="detect_free_surface_liquid_"]}]]
            set type [$node getAttribute "v"]
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="free_surface_factor_"]}]]
            set factor [$node getAttribute "v"]
            
            GiD_WriteCalculationFile puts {$$APPLY_DETECT_LIQUID_SURFACE} 
            if {$type == "No"} {
                GiD_WriteCalculationFile puts [= "0 %s" $factor]                                                
            } elseif {$type == "Yes"} {
                GiD_WriteCalculationFile puts [= "1 %s" $factor]
            }
        }
        
        
        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_material_model_solid_"]}]]
            set typemodel [$node getAttribute "v"]
            if {$typemodel == "Rigid body"} {
                set model "rigid_body"
            } elseif {$typemodel == "Linear Elasticity"} {
                set model "linear_elasticity"
            } elseif {$typemodel == "Mohr-Coulomb"} {
                set model "mohr_coulomb"
            } elseif {$typemodel == "External Material Model"} {
                set model "external_soil_model"}
            GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_SOLID}    
            GiD_WriteCalculationFile puts $model
            if {$typemodel == "Rigid body"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="x-constr"]}]]
                set type [$node getAttribute "v"]
                if {$type == "Yes"} {
                    set flag 1
                } else {
                    set flag 0
                }                        
                GiD_WriteCalculationFile puts {$$constraint_XDISPLACEMENT} 
                GiD_WriteCalculationFile puts $flag        
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="y-constr"]}]]
                set type [$node getAttribute "v"]
                if {$type == "Yes"} {
                    set flag 1
                } else {
                    set flag 0
                }                        
                GiD_WriteCalculationFile puts {$$constraint_YDISPLACEMENT} 
                GiD_WriteCalculationFile puts $flag
                if {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="z-constr"]}]]
                    set type [$node getAttribute "v"]
                    if {$type == "Yes"} {
                        set flag 1
                    } else {
                        set flag 0
                    }                        
                    GiD_WriteCalculationFile puts {$$constraint_ZDISPLACEMENT} 
                    GiD_WriteCalculationFile puts $flag        
                }                                                        
            } elseif {$typemodel == "Linear Elasticity"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_Young_modulus_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$YOUNG_MODULUS} 
                GiD_WriteCalculationFile puts $type
				set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$POISSON_RATIO} 
                GiD_WriteCalculationFile puts $type
                if {$typename == "Saturated material-undrained effective stress"} {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="undr_poisson_ratio_"]}]]
                    set type [$node getAttribute "v"]
                    GiD_WriteCalculationFile puts {$$POISSON_RATIO_UNDRAINED} 
                    GiD_WriteCalculationFile puts $type
                }
            } elseif {$typemodel == "Mohr-Coulomb"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_Young_modulus_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$YOUNG_MODULUS} 
                GiD_WriteCalculationFile puts $type                
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$POISSON_RATIO} 
                GiD_WriteCalculationFile puts $type 
				if {$typename == "Saturated material-undrained effective stress"} {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="undr_poisson_ratio_"]}]]
                    set type [$node getAttribute "v"]
                    GiD_WriteCalculationFile puts {$$POISSON_RATIO_UNDRAINED} 
                    GiD_WriteCalculationFile puts $type
                }
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_friction_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$FRICTION_ANGLE} 
                GiD_WriteCalculationFile puts $type                
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_cohesion_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$COHESION} 
                GiD_WriteCalculationFile puts $type        
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_dilatancy_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$DILATANCY_ANGLE} 
                GiD_WriteCalculationFile puts $type                
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_tensile_strength_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$TENSILE_STRENGTH} 
                GiD_WriteCalculationFile puts $type        
            } elseif {$typemodel == "External Material Model"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_model_dll_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_NAME} 
                GiD_WriteCalculationFile puts $type                
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_model_dll_dim_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$UMAT_DIMENSION} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/container[@n="mat_params"]/container[@n="list_material_parameter_state_var"]}]]
                set number_mat_paramsNode [$gNode selectNodes {container[@n="_material_constitutive_model"]/container[@n="mat_params"]/value[@n="number_mat_params"]}] 
                set number_mat_params [$number_mat_paramsNode @v]
                set number_state_varNode [$gNode selectNodes {container[@n="_material_constitutive_model"]/container[@n="mat_params"]/value[@n="number_state_var"]}]                                                
                set number_state_var [$number_state_varNode @v]                
                if {$node != ""} {
                    set icount 1
                    foreach iNode [$node selectNodes {value}] {
                        set type [$iNode @type]
                        if {$type != "number_mat_params"} { continue }
                        set vv [$iNode @v]
                        set nn [$iNode @n]
                        set nn [string range $nn 0 end-1]
                        set nn [string toupper $nn]
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv
                        incr icount                       
                    }
                    set var [regsub {_[0-9]+$} $nn ""]
                    for { set j $icount } { $j <= 50 } { incr j } {
                        if { $j<= 9 } {
                            set jj "0$j" 
                        } else {
                            set jj "$j"
                        }  
                        set nn ${var}_${jj}
                        set vv 0.0
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv                        
                    }                   
                    set icount 1
                    foreach iNode [$node selectNodes {value}] {
                        set type [$iNode @type]
                        if {$type != "number_state_var"} { continue } 
                        set vv [$iNode @v]
                        set nn [$iNode @n]
                        set nn [string range $nn 0 end-1]
                        set nn [string toupper $nn]
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv
                        incr icount                              
                    }
                    set var [regsub {_[0-9]+$} $nn ""]
                    for { set j $icount } { $j <= 50 } { incr j } {
                        if { $j<= 9 } {
                            set jj "0$j" 
                        } else {
                            set jj "$j"
                        }  
                        set nn ${var}_${jj}
                        set vv 0.0
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv                        
                    }                                                                         
                } else {                        
                    set var "MATERIAL_PARAMETER_SOLID"
                    for { set j 1 } { $j <= 50 } { incr j } {
                        if { $j<= 9 } {
                            set jj "0$j" 
                        } else {
                            set jj "$j"
                        }  
                        set nn ${var}_${jj}
                        set vv 0.0
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv                        
                    }
                    set var "INITIAL_STATE_VARIABLE_SOLID"
                    for { set j 1 } { $j <= 50 } { incr j } {
                        if { $j<= 9 } {
                            set jj "0$j" 
                        } else {
                            set jj "$j"
                        }  
                        set nn ${var}_${jj}
                        set vv 0.0
                        GiD_WriteCalculationFile puts $nn
                        GiD_WriteCalculationFile puts $vv                        
                    }
                }
            }   
        }                
        if {$typename == "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_material_model_liquid_"]}]]
            set typemodel [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_LIQUID}    
            if {$typemodel == "Newtonian"} {
                GiD_WriteCalculationFile puts  "newtonian_liquid"                                                                           
            } elseif {$typemodel == "Bingham Fluid"} {
                GiD_WriteCalculationFile puts  "bingham_liquid"
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="bingham_yield_stress_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_YIELD_STRESS} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_young_modulus_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_YOUNG_MODULUS}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_poisson_ratio_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_POISSON_RATIO}
                GiD_WriteCalculationFile puts $type
            } elseif {$typemodel == "Frictional Fluid"} {
                GiD_WriteCalculationFile puts  "frictional_liquid"
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="fluid_friction_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_FRICTION_ANGLE} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_young_modulus_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_YOUNG_MODULUS} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_poisson_ratio_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_POISSON_RATIO} 
                GiD_WriteCalculationFile puts $type
            }                                                                  
        }
        
        if {$typename == "Unsaturated material-2-phase with suction effect" || $typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_retention_curve_"]}]]
            set type [$node getAttribute "v"]
            if {$type == "Linear"} {
                GiD_WriteCalculationFile puts {$$WATER_RETENTION_CURVE}
                GiD_WriteCalculationFile puts {linear}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_rc_linear"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$av} 
                GiD_WriteCalculationFile puts $type        
            } elseif {$type == "Van Genuchten"} {
                GiD_WriteCalculationFile puts {$$WATER_RETENTION_CURVE}
                GiD_WriteCalculationFile puts {van_genuchten}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_min_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smin} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_max_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smax} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_ref_press_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$P0} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_lambda_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$L} 
                GiD_WriteCalculationFile puts $type        
            }
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_hydraulic_cond_"]}]]
            set type [$node getAttribute "v"]
            if {$type == "Constant"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {constant}      
            } elseif {$type == "Hillel"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {hillel}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_r_exp"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$r} 
                GiD_WriteCalculationFile puts $type
            } elseif {$type == "Mualem"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {mualem}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_min_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smin} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_max_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smax} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_ref_press_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$P0} 
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_lambda_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$L} 
                GiD_WriteCalculationFile puts $type        
            }
        }
    }
    
    # Material ID (2D/3D)        
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
        set ElementList [GiD_Info Mesh Elements Triangle -sublist]
        set list_len [llength $ElementList]
        set material_ID_list [lrepeat $list_len 0]
        set damping_list [lrepeat $list_len 0.0]
        set material_point_list_s [lrepeat $list_len 0]
        if {$layer_type == "Double_point"} {
            set material_point_list_l [lrepeat $list_len 0]
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Double-point"]/group} $ov_type]
        } else {
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Single-point"]/group} $ov_type]
        }        
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
        set ElementList [GiD_Info Mesh Elements Tetrahedra -sublist]
        set list_len [llength $ElementList]
        set material_ID_list [lrepeat $list_len 0]
        set damping_list [lrepeat $list_len 0.0]
        set material_point_list_s [lrepeat $list_len 0]
        if {$layer_type == "Double_point"} {
            set material_point_list_l [lrepeat $list_len 0]
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Double-point"]/group} $ov_type]
        } else {
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Single-point"]/group} $ov_type]
        }        
    }        
    foreach gNode [$stageNode selectNodes $xp] {       
        
        set l_material [$gNode selectNodes {string(value[@n="material"]/@v)}]
        set list_group [$gNode @n]
        set MATERIAL_ID [find_material_id $l_material $stageNode]
        set materialpoints_s [$gNode selectNodes {string(value[@n="solid_MP_number"]/@v)}]
        if {$layer_type == "Double_point"} {
            set materialpoints_l [$gNode selectNodes {string(value[@n="liquid_MP_number"]/@v)}]
        }
        set mat_damp [$gNode selectNodes {string(value[@n="material_damping"]/@v)}]
        set elements_id [GiD_EntitiesGroups get $list_group elements]
        set num_elems [objarray length $elements_id]
        for {set i 0} {$i < $num_elems} {incr i} {
            set node_id [objarray get $elements_id $i]
            lset material_ID_list [expr $node_id -1] $MATERIAL_ID
            lset material_point_list_s [expr $node_id -1] $materialpoints_s
            if {$layer_type == "Double_point"} {
                lset material_point_list_l [expr $node_id -1] $materialpoints_l
            }
            lset damping_list [expr $node_id -1] $mat_damp
        }
        
    }
    
    GiD_WriteCalculationFile puts {$$STARTELMMAT}
    set len_mat_id [llength $material_ID_list]
    for {set i 0} {$i < $len_mat_id} {incr i} {
        set print [lindex $material_ID_list $i]
        GiD_WriteCalculationFile puts $print 
    }
    
    GiD_WriteCalculationFile puts {$$STARTDAMPING}
    set len_dam_id [llength $damping_list]
    for {set i 0} {$i < $len_dam_id} {incr i} {
        set print [lindex $damping_list $i]
        GiD_WriteCalculationFile puts $print 
    }
    
    GiD_WriteCalculationFile puts {$$START_NUMBER_OF_MATERIAL_POINTS}
    if {$layer_type == "Single_point"} {
        set len_mps_id [llength $material_point_list_s]
        for {set i 0} {$i < $len_mps_id} {incr i} {
            set print [lindex $material_point_list_s $i]
            GiD_WriteCalculationFile puts [= "%s 0" $print]
        }
    } else {
        set len_mps_id [llength $material_point_list_s]
        for {set i 0} {$i < $len_mps_id} {incr i} {
            set prints [lindex $material_point_list_s $i]
            set printl [lindex $material_point_list_l $i]
            GiD_WriteCalculationFile puts [= "%s %s" $prints $printl]
        }        
    }

    GiD_WriteCalculationFile puts {$$FINISH}    
    GiD_WriteCalculationFile end
}
