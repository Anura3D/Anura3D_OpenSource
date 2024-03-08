# Write velocity information to the GOM file

proc write_prescribed_nodal_velocity {flag geometry_type dimension_type root_dir } {
    # Purpose: Write the prescribed nodal velocities for points, surfaces, and lines

    # Set the subpath for the prescribed velocity surface
    set sub_path [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $geometry_type]

    # set ov_type "point"
    # set sub_path [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="Nodal_prescribed_velocity"]/group[@ov=%s]} $ov_type]

    # Init dictionary to hold prescribed velocity data
    set presribed_vel_data ""

    # Select each node and store the data
    foreach gNode [$root_dir selectNodes $sub_path] {
	if {$dimension_type == "2D:plane-strain" || $dimension_type == "2D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
	    set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
	    set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
	    set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
		dict set $presribed_vel_data [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
	    } elseif {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
	    set v1 [$gNode selectNodes {string(value[@n="X_velocity_node"]/@v)}]
	    set v2 [$gNode selectNodes {string(value[@n="Y_velocity_node"]/@v)}]
	    set v3 [$gNode selectNodes {string(value[@n="Z_velocity_node"]/@v)}]
	    set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
	    set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
	    set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
		dict set $presribed_vel_data [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
	    }
    }

    # If data is collected then write the flag and the velocity data
    if { [dict size $presribed_vel_data]} {
	# Get the number of prescribed velocities
	set num [ GiD_WriteCalculationFile nodes -count $presribed_vel_data ]
	
	# Write the flag
	GiD_WriteCalculationFile puts $flag 

	# Write the number of prescribed velocites
	GiD_WriteCalculationFile puts $num

	# Write the prescribed velocities
	GiD_WriteCalculationFile puts $presribed_vel_data
    }

}

proc write_prescribed_MP_velocity { flag geometry_type root_dir } {
    # Purpose: Write the prescribed MP velocity for points, surfaces, and lines
    # Decided not to combine MP and nodal velocity procedures becasue MP have different conditons

    set sub_path [format_xpath {container[@n="BC"]/container[@n="Prescribed_velocity"]/condition[@n="MP_prescribed_velocity"]/group[@ov=%s]} $geometry_type]    

    # Init dictionary to hold the prescribed velocity data
    set prescribed_vel_data ""

    if {$geometry_type == "surface"} {
        # Select each surface node and store the data
        foreach gNode [$root_dir selectNodes $sub_path] {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_MP"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_MP"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            dict set prescribed_vel_data [$gNode @n] "%d $v1 $v2 $v3 $v4\n"
        }
    } elseif {$geometry_type == "volume"} {
        # Select each volume node and store the data
        foreach gNode [$root_dir selectNodes $sub_path] {
            set v1 [$gNode selectNodes {string(value[@n="X_velocity_MP"]/@v)}]
            set v2 [$gNode selectNodes {string(value[@n="Y_velocity_MP"]/@v)}]
            set v3 [$gNode selectNodes {string(value[@n="Z_velocity_MP"]/@v)}]
            set v4 [$gNode selectNodes {string(value[@n="X_velocity_[m/s]"]/@v)}]
            set v5 [$gNode selectNodes {string(value[@n="Y_velocity_[m/s]"]/@v)}]
            set v6 [$gNode selectNodes {string(value[@n="Z_velocity_[m/s]"]/@v)}]
            dict set prescribed_vel_data [$gNode @n] "%d $v1 $v2 $v3 $v4 $v5 $v6\n"
        }
    }
    # TODO: Add error conditon

    ## Write prescribed velocity data
    if { [dict size $prescribed_vel_data] } {
        set num [GiD_WriteCalculationFile elements -count $prescribed_vel_data]
        GiD_WriteCalculationFile puts $flag
        GiD_WriteCalculationFile puts $num
        GiD_WriteCalculationFile elements $prescribed_vel_data
    }
}