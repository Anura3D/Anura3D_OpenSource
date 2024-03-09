#      GOM file
#      print data in the .dat calculation file (instead of a classic .bas template)

# Construct the relative path to files that help create the GOM file
set general_output_path [file join [file dirname [info script]] GOMScripts writing_general_data.tcl]

# get path for writing element data
set element_output_path [file join [file dirname [info script]] GOMScripts write_element.tcl]

# get path for writing fixities
set fixity_output_path [file join [file dirname [info script]] GOMScripts write_fixities.tcl]

# Get the path for writing prescribed velocities
set presribed_velocity_path [file join [file dirname [info script]] GOMScripts write_velocities.tcl]

# Get the path for writing the material data
set material_output_path [file join [file dirname [info script]] GOMScripts write_material_data.tcl]

# Get the path for writing IGA data
set iga_output_path [file join [file dirname [info script]] GOMScripts write_iga_info.tcl]

# Include procedures from other files
source $general_output_path
source $element_output_path
source $fixity_output_path
source $presribed_velocity_path
source $material_output_path
source $iga_output_path

proc Anura3D::WriteCalculationFile_GOM { filename } {

    # Define some useful arrays
    # Set the geometry types
    set geometry_types(0) "POINT"
    set geometry_types(1) "LINE"
    set geometry_types(2) "SURFACE"
    set geometry_types(3) "VOLUMES"

    # Define an array of keys to hold 2D and lower geometries
    set geometry_types_2D_subset {0 1 2}

    # Set the material types
    set material_types(0) "solid"
    set material_types(1) "liquid"
    set material_types(2) "gas"

    # Set the dimension types
    set dimension_type_arr(0) "2D:plane-strain"
    set dimension_type_arr(1) "2D:Axissymmetric"
    set dimension_type_arr(2) "3D"
    set dimension_type_arr(3) "3D:Axissymmetric"

    # Init the file so the data is written to the correct one
    GiD_WriteCalculationFile init $filename

    set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
    set exe_name [GiD_Info Project ProblemType]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data
    customlib::SetBaseRoot $root
    set current_xml_root $root

    GiD_WriteCalculationFile puts "### Anura3D_2023 ###"

    ## Get the dimension from Gid
    set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
    set dim_type [$current_xml_root selectNodes $dim_path]

    # Get the geometry type
    set model_geometry_path {string(//container[@n = "Units_Dimensions"]/value[@n="model_geometry_type"]/@v)}
    set model_geometry_type [$current_xml_root selectNodes $model_geometry_path]

    # Get the formulation/layer selection path
    set layer_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
    set layer_type [$current_xml_root selectNodes $layer_path]

    # DIMENSION
    # Format the dimension and write it to GOM
    write_dimension {$$DIMENSION} $dim_type

    # MODEL GEOMETRY
    #Write the model geometry type and flag
    write_string_flag {$$MODEL_GEOMETRY_TYPE} $model_geometry_type

    # ELEMENT TYPE
    set info_mesh [GiD_Mesh get element 1]
    set elem_type [lindex $info_mesh  1]
    set elem_num_nodes [lindex $info_mesh  2]

    # Format and write the element type
    write_element_type {$$ELEMENTTYPE} $elem_type $elem_num_nodes

    # FORMULATION
    write_formulation_info {$$FORMULATION} $layer_type

    # Check if the geometry type is iga, if yes write the iga data
    if {$model_geometry_type eq "IGA"} {

        # Loop over the surfaces in the Gid_Geometry surfaces
        foreach surface_id [ GiD_Geometry list surface 1:end ] {

            # Store the surface data
            set surface_data [ GiD_Geometry get surface $surface_id ]

            # Format and write the nurbs surface data
            format_write_nurbs_data $surface_id $surface_data
        }
        # Writ the data for the regular geometry
    } elseif {$model_geometry_type eq "Regular"} {
        # COUNTERS
        set num_nodes [GiD_Info Mesh NumNodes]
        set num_elements [GiD_Info Mesh NumElements]
        write_reg_geometry_counters {$$STARTCOUNTERS} $num_elements $num_nodes

        # NODAL COORDINATES
        set Nodes [GiD_Info Mesh nodes -sublist]
        write_reg_geometry_nodes {$$STARTNODES} $num_nodes $Nodes $dim_type

        # ELEMENT CONNECTIVITIES
        write_element_connectivity {$$STARTELEMCON} $elem_type $elem_num_nodes $num_elements

        # FIXITIES
        ## Surface conditions
        ### Solid fixities
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="solid_fixities"]/group[@ov=%s]} $ov_type]
        # Set  varaible format to be an empty string
        set formats ""

        # Format and store solid surface fixities
        foreach gNode [$root selectNodes $xp] {
            if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
                set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
                set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
                dict set formats [$gNode @n] "%d $v1 $v2\n"
            } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
                set v1 [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
                set v2 [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
                set v3 [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write solid surface fixities to .GOM file
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid surface  fixities
        set ov_type "surface"
        set xp [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="liquid_fixities"]/group[@ov=%s]} $ov_type]
        set formats ""

        # Format and store liquid surface fixities
        foreach gNode [$root selectNodes $xp] {
            if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
                set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
                set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
                dict set formats [$gNode @n] "%d $v1 $v2\n"
            } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
                set v1 [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
                set v2 [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
                set v3 [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write liquid surface fixities to .GOM file
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_LIQUID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Gas surface fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write gas surface fixities to .GOM file
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_SURFACE_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ## Line Fixity conditions
        ### Solid line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write solid line fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write liquid line fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_LIQUID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Gas line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write Gas line fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_LINE_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ## Point conditions (Material point fixities)
        ### Solid MP fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ##### Write solid MP fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid MP fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write liquid MP fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_LIQUID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Gas MP fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write gas MP fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_FIXITY_POINT_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        # REMOVE FIXITIES
        ## Surface conditions

        ### Solid remove surface fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write Solid remove surface fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_SURFACE_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid remove surface fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write liquid remove surface fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write gas remove surface fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_SURFACE_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ## Line condtions removal
        ### Solid remove line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write solid remove line fixties
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid remove line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write liquid remove line fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_LIQUID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Gas remove line fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write gas remove line fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_LINE_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ## Point conditions
        ### Solid remove point fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write solid remove point fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_SOLID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Liquid remove point fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write liquid remove point fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_LIQUID}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats

        ### Gas remove point fixities
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        #### Write gas remove point fixities
        set num [GiD_WriteCalculationFile nodes -count $formats]
        GiD_WriteCalculationFile puts {$$START_REMOVE_FIXITY_POINT_GAS}
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile nodes $formats


        # NODAL PRESCRIBED VELOCITY (2D/3D)

        ## nodal prescribed velocity On point
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
            }
        }

        ### Write nodal prescribed velocity on point if there are values ($formats set in above foreach statement)
        if { [dict size $formats] } {
            set num [GiD_WriteCalculationFile nodes -count $formats]
            GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_POINT}
            GiD_WriteCalculationFile puts $num
            GiD_WriteCalculationFile nodes $formats
        }

        ## nodal prescribed velocity On line
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
            }
        }

        ### Write nodal prescribed velocity on line if there are values ($formats set in above foreach statement)
        if { [dict size $formats] } {
            et num [GiD_WriteCalculationFile nodes -count $formats]
            iD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_LINE}
            iD_WriteCalculationFile puts $num
            iD_WriteCalculationFile nodes $formats
        }

        ## nodal prescribed velocity On surface
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
            }
        }

        ### Write nodal prescribed velocity on surface if there are values ($formats set in above foreach statement)
        if { [dict size $formats] } {
            set num [GiD_WriteCalculationFile nodes -count $formats]
            GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_SURFACE}
            GiD_WriteCalculationFile puts $num
            GiD_WriteCalculationFile nodes $formats
        }

        ## nodal prescribed velocity On volume
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

        ### Write nodal prescribed velocity on volume if there are values ($formats set in above foreach statement)
        if { [dict size $formats] } {
            set num [GiD_WriteCalculationFile nodes -count $formats]
            GiD_WriteCalculationFile puts {$$PRESCRIBED_NODAL_VELOCITY_VOLUME}
            GiD_WriteCalculationFile puts $num
            GiD_WriteCalculationFile nodes $formats
        }

        # MATERIAL POINTS PRESCRIBED VELOCITY

        ## 2D - MP precribed velocity On surface
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

        ### Write 2D - MP prescribed velocity on surface
        if { [dict size $formats] } {
            set num [GiD_WriteCalculationFile elements -count $formats]
            GiD_WriteCalculationFile puts {$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE}
            GiD_WriteCalculationFile puts $num
            GiD_WriteCalculationFile elements $formats
        }

        # 3D - MP precribed velocity On volume
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

        ### Write 3D - MP precribed velocity on volume
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

        # TODO: Look into why select nodes is used as the iterating variables, it seems to be confusing
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3\n"
            }
        }

        ## Write initial velocity on material 2D/3D - (Writes material point initial velocities)
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
                set num_tot [expr $num + $num_tot]
            }
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
            dict set formats [$gNode @n] "%d $v1 \"$v2\" $v3 $v4 \"$v5\" $v6 $v7 \"$v8\" $v9 $v10 \"$v11\" $v12 $v13\n"
        }

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
            set excavation 1
        }

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
                    GiD_WriteCalculationFile puts [= "%s %s" $id_entity $id_element]
                }
            }
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
                GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax]
            }
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
                GiD_WriteCalculationFile puts [= "%s %s" $zmin $zmax]
            }
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
                GiD_WriteCalculationFile puts [= "%s %s %s" $x_value $y_value $z_value]
            }
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
        ## Surface conditions
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

        ## line absorbing boundary conditions
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
            }
        }
        set num2 [GiD_WriteCalculationFile nodes -count $formats]
        if {$num2 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_LINE_SOLID}
            GiD_WriteCalculationFile puts $num2
            GiD_WriteCalculationFile nodes $formats
        }
        ## node absorbing boundary conditions
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
            }
        }

        set num3 [GiD_WriteCalculationFile nodes -count $formats]

        if {$num3 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_POINT_SOLID}
            GiD_WriteCalculationFile puts $num3
            GiD_WriteCalculationFile nodes $formats
        }

        ## Surface absorbing boundary conditions
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
        # line abosorbing boundary conditions
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
            }
        }

        set num5 [GiD_WriteCalculationFile nodes -count $formats]
        if {$num5 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_LINE_LIQUID}
            GiD_WriteCalculationFile puts $num5
            GiD_WriteCalculationFile nodes $formats
        }

        ## node abosorbing boundary conditions
        ### Solid - node abosorbing boundary
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
                dict set formats [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9\n"
            }
        }

        set num6 [GiD_WriteCalculationFile nodes -count $formats]

        if {$num6 != 0} {
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_POINT_LIQUID}
            GiD_WriteCalculationFile puts $num6
            GiD_WriteCalculationFile nodes $formats
        }
        if {$num1 !=0 || $num2 !=0 || $num3 !=0 || $num4 !=0 || $num5 !=0 || $num6 !=0} {
            set xp [format_xpath {container[@n="BC"]/container[@n="Absorbing_boundary"]/value[@n="material_absorbing"] } ]
            set node [$root selectNodes $xp]
            set material_name [$node getAttribute "v"]
            set MATERIAL_ID [find_material_id $material_name $root]
            GiD_WriteCalculationFile puts {$$ABSORBING_BOUNDARY_REFERENCE_MATERIAL_INDEX}
            GiD_WriteCalculationFile puts $MATERIAL_ID
        }

        # Loading_Conditions (2D/3D)
        # 2D On line
        if { ($dim_type == "2D:plane-strain") || ($dim_type == "2D:Axissymmetric") } {
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
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
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
        } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
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
                        set flag 1
                    }
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
                        set flag 1
                    }
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
        ## Extending_mesh
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

        ## Compressing_mesh
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

        ## Moving_mesh
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

            ### Moving_mesh Reference material
            set xp [format_xpath {container[@n="Moving_mesh"]/value[@n="Reference_material"]}]
            set node [$root selectNodes $xp]
            set material_name [$node getAttribute "v"]
            set MATERIAL_ID [find_material_id $material_name $root]
            GiD_WriteCalculationFile puts {$$MOVING_MESH_REFERENCE_MATERIAL_INDEX}
            GiD_WriteCalculationFile puts $MATERIAL_ID
        }
    } else {
        error ["Geometry must be Regular or IGA"]
        # End condtion that checks if the geometry is regular or IGA and doesn't use the regular geometry boundary conditions modules
    }

    # MATERIALS
    # Write the material properties
    write_material_properties $dim_type $root

    # Material ID (2D/3D)
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
        if {$elem_type == "Triangle"} {
            set ElementList [GiD_Info Mesh Elements Triangle -sublist]
        } elseif {$elem_type == "Quadrilateral"} {
            set ElementList [GiD_Info Mesh Elements Quadrilateral -sublist]
        }
        set list_len [llength $ElementList]

        ## init material ids list with value of zero and length $list_len
        set material_ID_list [lrepeat $list_len 0]

        # Init damping factor list with value of 0.0
        set damping_list [lrepeat $list_len 0.0]

        # Init a list of the number of SOLID material points with value 0
        set material_point_list_s [lrepeat $list_len 0]

        # Set @n parameter and init list to hold liquid material points if double point formulation
        if {$layer_type == "Double_point"} {
            set material_point_list_l [lrepeat $list_len 0]
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Double-point"]/group} $ov_type]
        } else {
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Single-point"]/group} $ov_type]
        }
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
        if {$elem_type == "Tetrahedra"} {
            set ElementList [GiD_Info Mesh Elements Tetrahedra -sublist]
        } elseif {$elem_type == "Hexahedral"} {
            set ElementList [GiD_Info Mesh Elements Hexahedral -sublist]
        }

        set list_len [llength $ElementList]

        # init material ids list with value of zero and length $list_len
        set material_ID_list [lrepeat $list_len 0]

        # Init damping factor list with value of 0.0
        set damping_list [lrepeat $list_len 0.0]

        # Init a list of the number of SOLID material points with value 0
        set material_point_list_s [lrepeat $list_len 0]

        # Set @n parameter and init list to hold liquid material points if double point formulation
        if {$layer_type == "Double_point"} {
            set material_point_list_l [lrepeat $list_len 0]
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Double-point"]/group} $ov_type]
        } else {
            set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Single-point"]/group} $ov_type]
        }
    }

    foreach gNode [$root selectNodes $xp] {
        # List of materials (eg. SOIL_1)
        set l_material [$gNode selectNodes {string(value[@n="material"]/@v)}]

        # Grabs a material point specification object from the varaible $xp
        set list_group [$gNode @n]

        # material id associated with
        set MATERIAL_ID [find_material_id $l_material $root]

        set materialpoints_s [$gNode selectNodes {string(value[@n="solid_MP_number"]/@v)}]
        if {$layer_type == "Double_point"} {
            set materialpoints_l [$gNode selectNodes {string(value[@n="liquid_MP_number"]/@v)}]
        }

        set mat_damp [$gNode selectNodes {string(value[@n="material_damping"]/@v)}]

        # FIXME: The problem is $list_group elements isn't returning any element ids
        set elements_id [GiD_EntitiesGroups get $list_group elements]

        # FIXME: the length of $elements_id is zero causing num_elems (number of elements in a )
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

