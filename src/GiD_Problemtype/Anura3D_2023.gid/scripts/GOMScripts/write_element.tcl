proc write_element_type { element_flag element_geometry num_nodes } {
    # Purpose: Write the element type to the GOM file

    # Write the element type flag
    GiD_WriteCalculationFile puts $element_flag

    # Set the element based on the element geometry and the number of nodes
    if { ($element_geometry == "Triangle") && ($num_nodes == "3") } {
        GiD_WriteCalculationFile puts "triangular_3-noded"
    } elseif { ($element_geometry == "Quadrilateral") && ($num_nodes == "4") } {
        GiD_WriteCalculationFile puts "quadrilateral_4-noded"
    } elseif { ($element_geometry == "Quadrilateral") && ($num_nodes == "8") } {
        GiD_WriteCalculationFile puts "quadrilateral_8-noded"
    } elseif { ($element_geometry == "Tetrahedra") && ($num_nodes == "10") } {
        GiD_WriteCalculationFile puts "tetrahedral_old"
    } else {
        error [= "INPUT ERROR: Element type not properly defined. Only the following element types are supported: triangular 3-noded,
                  tetrahedral 10-noded, quadrilateral 4-noded, quadrilateral 8-noded."]
    }
}

proc write_triangle_connectivity { num_element_nodes num_elements element_list } {
    # Purpose: Write the connectivity for triangular elements

    if {$num_element_nodes == "3"} {
       for {set i 0} {$i < $num_elements } {incr i} {
           set xcoor [lindex $element_list  $i 1]
           set ycoor [lindex $element_list  $i 2]
           set zcoor [lindex $element_list  $i 3]
           GiD_WriteCalculationFile puts [= "%s %s %s" $xcoor $ycoor $zcoor]
       }
    } else {
        error [= "INPUT ERROR: Only 3 noded triangles are allowed"]
    }
}

proc write_quad_connectivity { num_element_nodes num_elements element_list } {
    # Purpose: Write the connectivity for quad elements

    # Select number of element nodes for quad elements and print to .GOM file
    if {$num_element_nodes == "4"} {
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $element_list  $i 1]
            set ycoor [lindex $element_list  $i 2]
            set zcoor [lindex $element_list  $i 3]
            set wcoor [lindex $element_list  $i 4]
            GiD_WriteCalculationFile puts [= "%s %s %s %s" $xcoor $ycoor $zcoor $wcoor]
        }
    } elseif {$num_element_nodes == "8"} {
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $element_list  $i 1]
            set ycoor [lindex $element_list  $i 2]
            set zcoor [lindex $element_list  $i 3]
            set wcoor [lindex $element_list  $i 4]
            set jcoor [lindex $element_list  $i 5]
            set qcoor [lindex $element_list  $i 6]
            set kcoor [lindex $element_list  $i 7]
            set lcoor [lindex $element_list  $i 8]
            GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s %s %s" $xcoor $ycoor $zcoor $wcoor $jcoor $qcoor $kcoor $lcoor]
        }
    }
}

proc write_tetrahedra_connectivity { num_element_nodes num_elements element_list } {
    # Purpose: Write the connectivity for tetrahedra elements
    
    # Select number of  element nodes for tetrahedra (3D triangle) and print to .GOM file
    if {$num_element_nodes == "10"} {
        for {set i 0} {$i < $num_elements } {incr i} {
            set xcoor [lindex $element_list  $i 1]
            set ycoor [lindex $element_list  $i 2]
            set zcoor [lindex $element_list  $i 3]
            set wcoor [lindex $element_list  $i 4]
            set jcoor [lindex $element_list  $i 5]
            set qcoor [lindex $element_list  $i 6]
            set kcoor [lindex $element_list  $i 7]
            set lcoor [lindex $element_list  $i 8]
            set mcoor [lindex $element_list  $i 9]
            set ncoor [lindex $element_list  $i 10]
            GiD_WriteCalculationFile puts [= "%s %s %s %s %s %s %s %s %s %s" $xcoor $ycoor $zcoor $wcoor $jcoor $qcoor $kcoor $lcoor $mcoor $ncoor]
        }
    }
}

proc write_element_connectivity { flag element_geometry num_element_nodes num_elements } {
    # Purpose: Act as a wrapper for the specific element print and write the flag

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Write the connectivities
    if {$element_geometry == "Triangle"} {
        # Get the element connectivity
        set element_list [GiD_Info Mesh Elements Triangle -sublist]

        # Write the triangle connectivity
        write_triangle_connectivity $num_element_nodes $num_elements $element_list

    } elseif {$element_geometry == "Quadrilateral"} {
        # Get the element connectivity
        set element_list [GiD_Info Mesh Elements Quadrilateral -sublist]

        # Write the quadrilateral connectivity
        write_quad_connectivity $num_element_nodes $num_elements $element_list

    } elseif {$element_geometry == "Tetrahedra"} {
        # Get the element connectivity
        set element_list [GiD_Info Mesh Elements Tetrahedra -sublist]

        # Write the tetrahedra connectivity
        write_tetrahedra_connectivity $num_element_nodes $num_elements $element_list
    } 
    # TODO: Add an error message here
}