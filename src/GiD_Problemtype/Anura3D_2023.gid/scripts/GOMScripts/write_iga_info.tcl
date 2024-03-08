proc write_multival_string {flag val_1 val_2} {
    # Purpose: Write val_1 and val_2 to the file

    # Write flag
    GiD_WriteCalculationFile puts $flag

    # Write the information
    GiD_WriteCalculationFile puts [concat $val_1 "" $val_2]
}

proc write_array_flag { flag array } {
    # Purpose: Write a an array under a flag
    
    # Write the flag
    GiD_WriteCalculationFile puts $flag

    foreach val $array {
	# Write each value on a new line
	GiD_WriteCalculationFile puts $val
    }
} 

# Function to flip the order of each array in the nested list
proc flipNestedList {nestedList} {
    set flippedList [list]
    foreach subList $nestedList {
	lappend flippedList [lreverse $subList]
    }
    return $flippedList
}

# Procedure to sort list in descending order and store mapping
proc sort_ListDescendingWithMapping {inputList} {
    set sortedList [lsort -decreasing $inputList]
    set mapping [list]

    foreach value $inputList {
        lappend mapping [lsearch -exact $sortedList $value]
    }

    return [list $sortedList $mapping]
}

# Procedure to reorder a list based on a mapping
proc reorderListUsingMapping {inputList mapping} {
    set reorderedList [list]

    foreach index $mapping {
        lappend reorderedList [lindex $inputList $index]
    }

    return $reorderedList
}

proc format_write_nurbs_data { surface_id surface_data } {
    # Pupose: Write and format the nurbs surface data
    # Check if the surface is a nurbs surface
    if { [lindex $surface_data 0] == "nurbssurface"} {
	
	# Get the innfo about the surface
	# Create an 8 element list and store the data
	lassign [lrange $surface_data 1 8] s_layer s_num_curves s_degree_u s_degree_v s_num_control_points_u s_num_control_points_v s_is_trimmed s_is_rational
	
	# Create a variables to track the number of curves
	set index_end [expr {9+$s_num_curves-1}]
	
	# Store the curve data
	set s_curves [lrange $surface_data 9 $index_end]
	
	# Store the geometry data
	set s_geometry [lrange $surface_data 9+$s_num_curves end]
	
	# Store the number of control points
	set s_num_control_points [expr {$s_num_control_points_u*$s_num_control_points_v}]
	
	# Set the index to be the number of control points my 1 (This will be the last element in the following arrays)
	set index_end [expr {$s_num_control_points-1}]
	
	# Get the control points data
	set s_control_points [lrange $s_geometry 0 $index_end]
	
	# Get the num of knots in both directions
	set s_num_knots_u [expr {$s_num_control_points_u+$s_degree_u+1}]
	set s_num_knots_v [expr {$s_num_control_points_v+$s_degree_v+1}]
	
	# Set the index to be the 
	set index_end [expr {$s_num_control_points+$s_num_knots_u-1}]
	set s_knots_u [lrange $s_geometry $s_num_control_points $index_end]
	set s_knots_v [lrange $s_geometry $index_end+1 end]  
	
	# Init list to store the control points in each direction                      
	set x_control_points [list]
	set y_control_points [list]
	set z_control_points [list]
	set control_point_weights [list]

	# Loop over the control points and store the directional data in arrays
	foreach control_point $s_control_points {
	    lassign $control_point x y z w
	    if { $w == "" } {
		set w 0 ;#if non rational w is missing, set to 0
	    }
	    lappend x_control_points $x
	    lappend y_control_points $y
	    lappend z_control_points $z                
	    lappend control_point_weights $w                
	}
	
	# Set variable to track to flip control points data
	set descending_order "1"

	if {$descending_order} {
        # Base the ordering on the x_control_points
        set result [sort_ListDescendingWithMapping $x_control_points ]

        set x_control_points [lindex $result 0]
        set mapping [lindex $result 1]

	    # # rearrange the order of the other lists using the mapping
	    set y_control_points [reorderListUsingMapping $y_control_points $mapping]
	    set z_control_points [reorderListUsingMapping $z_control_points $mapping]
	    set control_point_weights [reorderListUsingMapping $control_point_weights $mapping]
	}

	# Need the stuff below here
	# WRITING DATA

	# Write patch id
	write_string_flag {$$PATCH_ID} $surface_id

	# Write the number of knots and the order
	write_multival_string {$$XI_NUMBER_OF_KNOTS_AND_ORDER} $s_num_knots_u $s_degree_u
	write_multival_string {$$ETA_NUMBER_OF_KNOTS_AND_ORDER} $s_num_knots_v $s_degree_v

	# Write the number of control points
	write_string_flag {$$NUMBER_OF_CONTROL_POINTS} $s_num_control_points

	#Write the xi and eta knots
	write_array_flag {$$XI_KNOT} $s_knots_u
	write_array_flag {$$ETA_KNOT} $s_knots_v

    GiD_WriteCalculationFile puts {$$NURBS_CONTROL_POINTS}
    # Loop over the elements of the control points coordinates list
    foreach index [range 0 [expr {[llength $x_control_points] - 1}]] {
        set x_val [lindex $x_control_points $index]
        set y_val [lindex $y_control_points $index]
        set z_val [lindex $z_control_points $index]

        # Write the values to the GOM file
        GiD_WriteCalculationFile puts [concat $x_val $y_val $z_val]
    }


	# # Write the weights
	write_array_flag {$$NURBS_CONTROL_POINT_WEIGHTS } $control_point_weights

	# Write the number of knots and order of traction
    write_multival_string {$$XI_NUMBER_OF_KNOTS_AND_ORDER_TRACTION} 0 0
    write_multival_string {$$ETA_NUMBER_OF_KNOTS_AND_ORDER_TRACTION} 0 0

	# Number of traction control points
    write_string_flag {$$NUMBER_OF_CONTROL_POINTS_TRACTION} 0
	
    # xi and eta traction knots
    write_string_flag {$$XI_KNOT_TRACTION} 0

    write_string_flag {$$ETA_KNOT_TRACTION} 0

	# Write the traction control points
    write_string_flag {$$NURBS_CONTROL_POINTS_TRACTION} 0

	# Write the control points weights
    write_string_flag {$$NURBS_CONTROL_POINT_WEIGHTS_TRACTION} 0

	# Write which local elements the traction is applied to 
    write_string_flag {$$TRACTION_APPLY_LOCAL_ELEMENTS} 0
    }
}


