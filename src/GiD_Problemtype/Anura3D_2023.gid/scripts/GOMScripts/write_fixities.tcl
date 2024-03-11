# Procedures for formatting and writting fixity data
# TODO: Collpase the functions for the different materials into the same one (maybe? this also might be ok)

proc write_solid_fixity { flag geometry_type dimension_type root_dir } {
    # Purpose write solid fixities

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Start to collect the fixity data...

    # Get the subpath to the fixity data
    set sub_path [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="solid_fixities"]/group[@ov=%s]} $geometry_type]

    # Inot dictionary to hold fixitity data
    set fixity_data ""

    # Format and store solid "geometry_type" fixities
    foreach gNode [$root_dir selectNodes $sub_path] {
        if {$dimension_type == "2D:plane-strain" || $dimension_type == "2D:Axissymmetric"} {

            # Store the x and y constraints
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]

            # Store constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint\n"
            #
        } elseif {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint"]/@v)}]
            set z_constraint [$gNode selectNodes {string(value[@n="Z_Constraint"]/@v)}]

            # Store the constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint $z_constraint\n"
        }
    }

    # Get the number of fixities
    set num_fixity [GiD_WriteCalculationFile nodes -count $fixity_data]

    # Write number of and the actual solid "geometry_type" fixities
    GiD_WriteCalculationFile puts $num_fixity
    GiD_WriteCalculationFile nodes $fixity_data
}

proc write_liquid_fixity { flag geometry_type dimension_type root_dir } {
    # Purpose write liquid fixities

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Start to collect the fixity data...

    # Get the subpath to the fixity data
    set sub_path [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="liquid_fixities"]/group[@ov=%s]} $geometry_type]

    # Inot dictionary to hold fixitity data
    set fixity_data ""

    # Format and store liquid "geometry_type" fixities
    foreach gNode [$root_dir selectNodes $sub_path] {
        if {$dimension_type == "2D:plane-strain" || $dimension_type == "2D:Axissymmetric"} {

            # Store the x and y constraints
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]

            # Store constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint\n"
            #
        } elseif {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint_liq"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint_liq"]/@v)}]
            set z_constraint [$gNode selectNodes {string(value[@n="Z_Constraint_liq"]/@v)}]

            # Store the constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint $z_constraint\n"
        }
    }

    # Get the number of fixities
    set num_fixity [GiD_WriteCalculationFile nodes -count $fixity_data]

    # Write number of and the actual liquid "geometry_type" fixities
    GiD_WriteCalculationFile puts $num_fixity
    GiD_WriteCalculationFile nodes $fixity_data
}

proc write_gas_fixity { flag geometry_type dimension_type root_dir } {
    # Purpose write gas fixities

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Start to collect the fixity data...

    # Get the subpath to the fixity data
    set sub_path [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n="gas_fixities"]/group[@ov=%s]} $geometry_type]

    # Inot dictionary to hold fixitity data
    set fixity_data ""

    # Format and store gas "geometry_type" fixities
    foreach gNode [$root_dir selectNodes $sub_path] {
        if {$dimension_type == "2D:plane-strain" || $dimension_type == "2D:Axissymmetric"} {

            # Store the x and y constraints
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]

            # Store constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint\n"
            #
        } elseif {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
            set x_constraint [$gNode selectNodes {string(value[@n="X_Constraint_gas"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_Constraint_gas"]/@v)}]
            set z_constraint [$gNode selectNodes {string(value[@n="Z_Constraint_gas"]/@v)}]

            # Store the constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint $z_constraint\n"
        }
    }

    # Get the number of fixities
    set num_fixity [GiD_WriteCalculationFile nodes -count $fixity_data]

    # Write number of and the actual gas "geometry_type" fixities
    GiD_WriteCalculationFile puts $num_fixity
    GiD_WriteCalculationFile nodes $fixity_data
}

proc write_remove_fixity {flag geometry_type dimension_type phase root_dir} {
    # Purpose: Write the remove fixities conditions for solid, liquid, and gas

    # Write the flag
    GiD_WriteCalculationFile puts $flag

    # Get the phase string for selecting the correct fixities
    if { $phase eq "solid" } {
        set condition_phase_tag "solid"
        set constraint_phase_tag ""

    } elseif { $phase eq "liquid" } {
        set condition_phase_tag "liquid"
        set constraint_phase_tag "_liq"
    } elseif { $phase eq "gas" } {
        set condition_phase_tag "gas"
        set constraint_phase_tag "_gas"
    } else {
        error ["Improper phase selection for remove fixities. Should be solid, liquid, or gas"]
    }

    # Construct the condtion and the contraint values
    set condition [concat "remove_$condition_phase_tag\_fixities"]
    set constraint [concat "Constraint$constraint_phase_tag"]

    # Get the subpath to the fixity data
    set sub_path [format_xpath {container[@n="BC"]/container[@n="Fixities"]/condition[@n=$condition]/group[@ov=%s]} $geometry_type]

    # Inot dictionary to hold fixitity data
    set fixity_data ""

    # Select the fixities
    foreach gNode [$root_dir selectNodes $sub_path] {
        if {$dimension_type == "2D:plane-strain" || $dimension_type == "2D:Axissymmetric"} {

            # Store the x and y constraints
            set x_constraint [$gNode selectNodes {string(value[@n="X_$constraint"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_$constraint"]/@v)}]

            # Store constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint\n"
            #
        } elseif {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
            set x_constraint [$gNode selectNodes {string(value[@n="X_$constraint"]/@v)}]
            set y_constraint [$gNode selectNodes {string(value[@n="Y_$constraint"]/@v)}]
            set z_constraint [$gNode selectNodes {string(value[@n="Z_$constraint"]/@v)}]

            # Store the constraints in dictionary
            dict set fixity_data [$gNode @n] "%d $x_constraint $y_constraint $z_constraint\n"
        }
    }

    # Get the number of fixities
    set num_fixity [GiD_WriteCalculationFile nodes -count $fixity_data]

    # Write number of and the actual gas "geometry_type" fixities
    GiD_WriteCalculationFile puts $num_fixity
    GiD_WriteCalculationFile nodes $fixity_data

}