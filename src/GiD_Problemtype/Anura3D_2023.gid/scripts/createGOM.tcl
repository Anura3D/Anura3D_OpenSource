#      GOM file
#      print data in the .dat calculation file (instead of a classic .bas template)


proc Anura3D::WriteCalculationFile_GOM { filename } {
  variable current_xml_root
  GiD_WriteCalculationFile init $filename
    set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
    set exe_name [GiD_Info Project ProblemType]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data
    customlib::SetBaseRoot $root
    set current_xml_root $root

  GiD_WriteCalculationFile puts "### Anura3D_2023 ###"
    
      # DIMENSION
      GiD_WriteCalculationFile puts {$$DIMENSION}
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
      set dim_type [$current_xml_root selectNodes $dim_path]
	if {$dim_type == "2D:plane-strain"} {
	  GiD_WriteCalculationFile puts "2D-plane_strain"
	} elseif {$dim_type == "2D:Axissymmetric"} {
	  GiD_WriteCalculationFile puts "2D-axisymmetric"
	} elseif {$dim_type == "3D"} {
	  GiD_WriteCalculationFile puts "3D-cartesian"
	} elseif {$dim_type == "3D:Axissymmetric"} {
	  GiD_WriteCalculationFile puts "3D-cylindrical"}

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
      set layer_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
      set layer_type [$current_xml_root selectNodes $layer_path]
	if {$layer_type == "Single_point"} {
	  GiD_WriteCalculationFile puts "single-point"
	} elseif {$layer_type == "Double_point"} {
	  GiD_WriteCalculationFile puts "double-point"}

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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {      
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
	foreach gNode [$root selectNodes $xp] {      
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
	foreach gNode [$root selectNodes $xp] {      
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
    
	# REMOVE FIXITIES
	# Surface conditions
	### Solid
	set ov_type "surface"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_solid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
		if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
		if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
		if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
		if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
		if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_SURFACE_SOLID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
	### Liquid
	set ov_type "surface"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_liquid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_SURFACE_LIQUID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
	### Gas
	set ov_type "surface"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_gas_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_SURFACE_GAS}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats 
    
	# Line conditions
	### Solid
	set ov_type "line"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_solid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {      
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}        
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_SOLID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
	### Liquid
	set ov_type "line"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_liquid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {      
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}        
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_LIQUID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats 
	### Gas
	set ov_type "line"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_gas_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {      
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}        
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_GAS}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
	
	# Point conditions
	### Solid
	set ov_type "point"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_solid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {       
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_SOLID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
	### Liquid
	set ov_type "point"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_liquid_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {       
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_LIQUID}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats       
	### Gas
	set ov_type "point"
	set xp [format_xpath {container[@n="BC"]/container[@n="Remove_Fixities"]/condition[@n="remove_gas_fixities"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {       
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    dict set formats [$gNode @n] "%d $v1 $v2\n"
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
	    if {$v1 == "1"} {set v1 -1}
	    set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
	    if {$v2 == "1"} {set v2 -1}
	    set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]
	    if {$v3 == "1"} {set v3 -1}
	     dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"}
	}
	set num [GiD_WriteCalculationFile nodes -count $formats]
	GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_GAS}    
	GiD_WriteCalculationFile puts $num
	GiD_WriteCalculationFile nodes $formats
    
    
	# NODAL PRESCRIBED VELOCITY (2D/3D)
	# On point
	set ov_type "point"
	set xp [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {       
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
	foreach gNode [$root selectNodes $xp] {
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
	set xp [format_xpath {container[@n="Initial_cond"]/container[@n="Phreatic_surface"]/blockdata}]
	set list [$root selectNodes $xp]        
	set num_tot 0
    
	foreach gNode $list { 
	set num [$gNode selectNodes [format_xpath {string(value[@n="Number_of_materials"]/@v)}]]
	if {$num == "please specify"} {
	set num 0
	} else {
	set num_tot [expr $num + $num_tot] }
	}
    
	if {$num_tot != 0} {
	    GiD_WriteCalculationFile puts {$$INITIAL_WATER_SURFACE_FROM_FILE}
	    GiD_WriteCalculationFile puts $num_tot
	}
    
	foreach gNode $list {
	set name [$gNode getAttribute "name"]
	set xp [format_xpath {string(value[@n="Phreatic_surface_file_flag"]/@v)}]
	set type_flag [$gNode selectNodes $xp]
	if {$type_flag == "yes"} {
	    if {$name == "Water table from file PSF_1"} {
	    set type_number 1
	    } elseif {$name == "Water table from file PSF_2"} {        
	    set type_number 2
	    } elseif {$name == "Water table from file PSF_3"} {
	    set type_number 3
	}      
	set num [$gNode selectNodes [format_xpath {string(value[@n="Number_of_materials"]/@v)}]]     
	set type_name [$gNode selectNodes [format_xpath {string(value[@n="material_phreatic_surface_file1"]/@v)}]]      
	set MATERIAL_ID [find_material_id $type_name $root]
	    GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
	if {$num >= "2"} {
	    set type_name [$gNode selectNodes [format_xpath {string(value[@n="material_phreatic_surface_file2"]/@v)}]]        
	    set MATERIAL_ID [find_material_id $type_name $root]
	    GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
	}        
	if {$num >= "3"} {
	    set type_name [$gNode selectNodes [format_xpath {string(value[@n="material_phreatic_surface_file3"]/@v)}]]
	    set MATERIAL_ID [find_material_id $type_name $root]
	    GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
	} 
	if {$num >= "4"} {
	    set type_name [$gNode selectNodes [format_xpath {string(value[@n="material_phreatic_surface_file4"]/@v)}]]
	    set MATERIAL_ID [find_material_id $type_name $root]
		GiD_WriteCalculationFile puts [= "%s %s" $MATERIAL_ID $type_number]
	} 
	if {$num == "5"} { 
	    set type_name [$gNode selectNodes [format_xpath {string(value[@n="material_phreatic_surface_file5"]/@v)}]]
	    set MATERIAL_ID [find_material_id $type_name $root]
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
	set xp [format_xpath {condition[@n="Contact_properties"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {       
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
    foreach gNode [$root selectNodes $xp] {
	set excavation 1}
	if {$excavation == 1} {
	GiD_WriteCalculationFile puts {$$START_EXCAVATION_SOLID}
	}
	foreach gNode [$root selectNodes $xp] {
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
	set hydraulic_head 0
	set ov_type "point"        
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/condition[@n="Minimum_coordinates_hydraulic_head_2D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/condition[@n="Maximum_coordinates_hydraulic_head_2D"]/group[@ov=%s]} $ov_type]
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/condition[@n="Minimum_coordinates_hydraulic_head_3D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/condition[@n="Maximum_coordinates_hydraulic_head_3D"]/group[@ov=%s]} $ov_type]    
	}
	foreach gNode [$root selectNodes $xp_min] {
	    set hydraulic_head 1
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmin [lindex $node_coord  0]
	    set ymin [lindex $node_coord  1] 
	    set zmin [lindex $node_coord  2] 
	}
	foreach gNode [$root selectNodes $xp_max] {
	    set hydraulic_head 2
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmax [lindex $node_coord  0]
	    set ymax [lindex $node_coord  1]
	    set zmax [lindex $node_coord  2] 
	}
	if {$hydraulic_head == 2} {
	    GiD_WriteCalculationFile puts {$$BOUNDARY_HYDRAULIC_HEAD_AREA}
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
	    GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {  
	    GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
	    GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
	    GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax] }
	}
	 
	# Seepage face
	set seepage_face 0
	set ov_type "point"
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/condition[@n="Minimum_coordinates_seepage_face_2D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/condition[@n="Maximum_coordinates_seepage_face_2D"]/group[@ov=%s]} $ov_type]
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/condition[@n="Minimum_coordinates_seepage_face_3D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/condition[@n="Maximum_coordinates_seepage_face_3D"]/group[@ov=%s]} $ov_type]    
	}
	foreach gNode [$root selectNodes $xp_min] {
	    set seepage_face 1
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmin [lindex $node_coord  0]
	    set ymin [lindex $node_coord  1] 
	    set zmin [lindex $node_coord  2] 
	}
	foreach gNode [$root selectNodes $xp_max] {
	    set seepage_face 2
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmax [lindex $node_coord  0]
	    set ymax [lindex $node_coord  1]
	    set zmax [lindex $node_coord  2] 
	}
	if {$seepage_face == 2} {
	    GiD_WriteCalculationFile puts {$$BOUNDARY_SEEPAGE_AREA}
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
	    GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {  
	    GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
	    GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
	    GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax] }
	}    
	
	# Infiltration
	set infiltration 0
	set ov_type "point"
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/condition[@n="Minimum_coordinates_infiltration_2D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/condition[@n="Maximum_coordinates_infiltration_2D"]/group[@ov=%s]} $ov_type]
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	    set xp_min [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/condition[@n="Minimum_coordinates_infiltration_3D"]/group[@ov=%s]} $ov_type]
	    set xp_max [format_xpath {container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/condition[@n="Maximum_coordinates_infiltration_3D"]/group[@ov=%s]} $ov_type]    
	}
	foreach gNode [$root selectNodes $xp_min] {
	    set infiltration 1
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmin [lindex $node_coord  0]
	    set ymin [lindex $node_coord  1] 
	    set zmin [lindex $node_coord  2] 
	}
	foreach gNode [$root selectNodes $xp_max] {
	    set infiltration 2
	    set gName [get_domnode_attribute $gNode n]
	    set node_id [GiD_EntitiesGroups get $gName node]
	    set node_coord [GiD_Mesh get node $node_id coordinates]             
	    set xmax [lindex $node_coord  0]
	    set ymax [lindex $node_coord  1]
	    set zmax [lindex $node_coord  2] 
	}    
	if {$infiltration == 2} {    
	    set xp [format_xpath {string(//container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="X_direction"]/@v)}]
	    set x_value [$root selectNodes $xp]
	    set yp [format_xpath {string(//container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="Y_direction"]/@v)}]
	    set y_value [$root selectNodes $yp]
	    set zp [format_xpath {string(//container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/container[@n="Infiltration_rate_value"]/value[@n="Z_direction"]/@v)}]  
	    set z_value [$root selectNodes $zp]   
	    GiD_WriteCalculationFile puts {$$BOUNDARY_INFILTRATION_AREA}
	    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
		GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
		GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
		GiD_WriteCalculationFile puts {$$INFILTRATION_RATE}
		GiD_WriteCalculationFile puts [= "%s %s" $x_value $y_value]
	    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {  
		GiD_WriteCalculationFile puts [= "%s %s" $xmin $xmax]
		GiD_WriteCalculationFile puts [= "%s %s" $ymin $ymax]
		GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax] 
		GiD_WriteCalculationFile puts {$$INFILTRATION_RATE} 
		GiD_WriteCalculationFile puts [= "%s %s %s" $x_value $y_value $z_value] }
	}         
	
	# Reaction forces 2D/3D  (2D/3D)
	# 2D On line
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	set ov_type "line"
	set xp [format_xpath {condition[@n="Reaction_forces"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
		set line_name [$gNode selectNodes {string(value[@n="line_identifier"]/@v)}]
		dict set formats [$gNode @n] "\"$line_name\" %d %d %d\n"
	}    
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	set ov_type "surface"
	set xp [format_xpath {condition[@n="Reaction_forces"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
	
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	set node [$root selectNodes $xp]
	set material_name [$node getAttribute "v"]
	set MATERIAL_ID [find_material_id $material_name $root]
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
	foreach gNode [$root selectNodes $xp] {        

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
	foreach gNode [$root selectNodes $xp] {       

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
	foreach gNode [$root selectNodes $xp] {        

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
	foreach gNode [$root selectNodes $xp] {       

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
	set xp [format_xpath {container[@n="Initial_cond"]/condition[@n="Soil_surface"]/group[@ov=%s]} $ov_type]
	set formats ""       
	foreach gNode [$root selectNodes $xp] {
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
	set xp [format_xpath {container[@n="Initial_cond"]/container[@n="Phreatic_surface"]/condition[@n="Phreatic_surface_line"]/group[@ov=%s]} $ov_type]
	set formats ""
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	foreach gNode [$root selectNodes $xp] {
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
	set node [$root selectNodes $xp]
	set material_name [$node getAttribute "v"]
	set MATERIAL_ID [find_material_id $material_name $root]
	GiD_WriteCalculationFile puts {$$MOVING_MESH_REFERENCE_MATERIAL_INDEX}    
	GiD_WriteCalculationFile puts $MATERIAL_ID
	}

	# MATERIALS
	set xp [format_xpath {container[@n="materials"]/blockdata}]
	set list [$root selectNodes $xp]
	set list_len [llength $list]
		                GiD_WriteCalculationFile puts {$$NUMBER_OF_MATERIALS}    
		                GiD_WriteCalculationFile puts $list_len 
	set int 1
	foreach gNode $list {  
		
	set type [$gNode getAttribute "name"]
	set MATERIAL_ID [find_material_id $type $root]
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
	 if {$type == "yes"} {
		        set flag 1
	 } else {
		        set flag 0
	 }                        
		                GiD_WriteCalculationFile puts {$$constraint_XDISPLACEMENT} 
		                                GiD_WriteCalculationFile puts $flag        
	    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="y-constr"]}]]
		set type [$node getAttribute "v"]
	 if {$type == "yes"} {
		        set flag 1
	 } else {
		        set flag 0
	 }                        
		                GiD_WriteCalculationFile puts {$$constraint_YDISPLACEMENT} 
		                                GiD_WriteCalculationFile puts $flag
	if {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
		set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="z-constr"]}]]
		set type [$node getAttribute "v"]
	 if {$type == "yes"} {
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
		if {$typename != "Saturated material-undrained total stress"} {
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$POISSON_RATIO} 
		                                GiD_WriteCalculationFile puts $type
		} else {
		set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$POISSON_RATIO} 
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
		                        GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_DLL} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_model_dll_dim_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$UMAT_DIMENSION} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_01_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_01} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_02_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_02} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_03_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_03} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_04_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_04} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_05_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_05} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_06_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_06} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_07_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_07} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_08_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_08} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_09_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_09} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_10_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_10} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_11_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_11} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_12_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_12} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_13_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_13} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_14_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_14} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_15_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_15} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_16_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_16} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_17_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_17} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_18_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_18} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_19_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_19} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_20_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_20} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_21_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_21} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_22_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_22} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_23_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_23} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_24_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_24} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_25_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_25} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_26_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_26} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_27_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_27} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_28_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_28} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_29_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_29} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_30_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_30} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_31_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_31} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_32_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_32} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_33_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_33} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_34_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_34} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_35_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_35} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_36_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_36} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_37_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_37} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_38_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_38} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_39_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_39} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_40_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_40} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_41_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_41} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_42_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_42} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_43_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_43} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_44_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_44} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_45_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_45} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_46_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_46} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_47_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_47} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_48_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_48} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_49_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_49} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_50_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_50} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_01_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_01} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_02_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_02} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_03_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_03} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_04_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_04} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_05_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_05} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_06_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_06} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_07_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_07} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_08_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_08} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_09_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_09} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_10_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_10} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_11_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_11} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_12_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_12} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_13_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_13} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_14_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_14} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_15_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_15} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_16_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_16} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_17_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_17} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_18_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_18} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_19_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_19} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_20_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_20} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_21_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_21} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_22_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_22} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_23_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_23} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_24_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_24} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_25_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_25} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_26_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_26} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_27_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_27} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_28_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_28} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_29_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_29} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_30_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_30} 
		                                GiD_WriteCalculationFile puts $type                        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_31_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_31} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_32_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_32} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_33_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_33} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_34_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_34} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_35_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_35} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_36_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_36} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_37_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_37} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_38_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_38} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_39_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_39} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_40_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_40} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_41_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_41} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_42_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_42} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_43_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_43} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_44_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_44} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_45_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_45} 
		                                GiD_WriteCalculationFile puts $type                
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_46_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_46} 
		                                GiD_WriteCalculationFile puts $type
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_47_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_47} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_48_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_48} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_49_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_49} 
		                                GiD_WriteCalculationFile puts $type        
		 set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_50_"]}]]
		set type [$node getAttribute "v"]
		                        GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_50} 
		                                GiD_WriteCalculationFile puts $type                                                        
	 }         
	 }
	 if {$typename == "Liquid"} {
	    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_material_model_liquid_"]}]]
		set typemodel [$node getAttribute "v"]
		                GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_LIQUID}    
		                GiD_WriteCalculationFile puts $typemodel
	 if {$typemodel == "Newtonian"} {
	    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="x-constr"]}]]
		set type [$node getAttribute "v"]
		                GiD_WriteCalculationFile puts {$$constraint_XDISPLACEMENT} 
		                                GiD_WriteCalculationFile puts $type        

	if {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
		set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="z-constr"]}]]
		set type [$node getAttribute "v"]        
		                GiD_WriteCalculationFile puts {$$constraint_ZDISPLACEMENT} 
		                                GiD_WriteCalculationFile puts $flag        
	}                                                        
	 } elseif {$typemodel == "Bingham Fluid"} {
	 } elseif {$typemodel == "Frictional Fluid"} {
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
	foreach gNode [$root selectNodes $xp] {       
		
	set l_material [$gNode selectNodes {string(value[@n="material"]/@v)}]
	set list_group [$gNode @n]
		set MATERIAL_ID [find_material_id $l_material $root]
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
