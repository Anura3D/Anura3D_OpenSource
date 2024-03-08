proc write_string_flag { flag string } {
    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Write the flag data
    GiD_WriteCalculationFile puts $string
}

proc write_dimension { dimension_flag dimension_type } {
    # Purpose: Write the dimension of the model

    # Write the dimension flag
    GiD_WriteCalculationFile puts $dimension_flag
    
    # Select the dimension and write it
    if {$dimension_type == "2D:plane-strain"} {
        GiD_WriteCalculationFile puts "2D-plane_strain"
    } elseif {$dimension_type == "2D:Axissymmetric"} {
       GiD_WriteCalculationFile puts "2D-axisymmetric"
    } elseif {$dimension_type == "3D"} {
       GiD_WriteCalculationFile puts {"3D-cartesian"}
    } elseif {$dimension_type == "3D:Axissymmetric"} {
       GiD_WriteCalculationFile puts "3D-cylindrical"
    }
}

proc write_formulation_info { flag type } {
    # Purpose: Writ the fomulation of the model (Single or double point)

    GiD_WriteCalculationFile puts $flag

    # Select single point or double point formulation
    if {$type == "Single_point"} {
        GiD_WriteCalculationFile puts "single-point"
    } elseif {$type == "Double_point"} {
        GiD_WriteCalculationFile puts "double-point"
    }
}

proc write_reg_geometry_counters { flag num_elements num_nodes } {
    # Purpose: Write the number of elements and num of nodes for the regular geometry to serve as counters in Anura3D

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Write counters
    GiD_WriteCalculationFile puts [= "%s %s" $num_elements $num_nodes]
}

proc write_reg_geometry_nodes { flag num_nodes nodes dimension_type } {
    # Purpose: Write the nodes for the regular geometry mode

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Write the nodes
    if {$dimension_type == "2D:plane-strain"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $nodes  $i 1]
            set ycoor [lindex $nodes  $i 2]
            GiD_WriteCalculationFile puts [= "%s %s" $xcoor $ycoor]
        }
    } elseif {$dimension_type == "2D:Axissymmetric"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $nodes  $i 1]
            set ycoor [lindex $nodes  $i 2]
            GiD_WriteCalculationFile puts [= "%s %s" $xcoor $ycoor]
        }
    } elseif {$dimension_type == "3D"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $nodes  $i 1]
            set ycoor [lindex $nodes  $i 2]
            set zcoor [lindex $nodes  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]
        }
    } elseif {$dimension_type == "3D:Axissymmetric"} {
        for {set i 0} {$i < $num_nodes } {incr i} {
            set xcoor [lindex $nodes  $i 1]
            set ycoor [lindex $nodes  $i 2]
            set zcoor [lindex $nodes  $i 3]
            GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]
        }
    }
}



