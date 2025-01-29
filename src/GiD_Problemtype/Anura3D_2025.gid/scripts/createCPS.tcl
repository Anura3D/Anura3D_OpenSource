#      CPS file
#      print data in the .dat calculation file (instead of a classic .bas template)


proc Anura3D::WriteCalculationFile_CPS { filename stageNode icount_stage total_time number_of_active_elements out_num_material_points  out_material_points } {
    variable current_xml_root
       
    GiD_WriteCalculationFile init $filename
    set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
    set exe_name [GiD_Info Project ProblemType]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data
    set current_xml_root $root

    GiD_WriteCalculationFile puts "### Anura3D_2025 ###"
    
    # STAGE NUMBER
    GiD_WriteCalculationFile puts {$$STAGE}
    GiD_WriteCalculationFile puts $icount_stage        

    # PROBLEM DIMENSIONS
    set xp_dim {container[@n="General_data"]/value[@n="NDIM"]}
    set dim_typeNode [$current_xml_root selectNodes $xp_dim]
    set dim_type [$dim_typeNode @v]

    # NUMBER OF LOADSTEP
    GiD_WriteCalculationFile puts {$$NUMBER_OF_LOADSTEPS}
    set num_step_path {string(container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="number_of_calculation_steps"]/@v)}
    set num_step [$stageNode selectNodes $num_step_path]
    GiD_WriteCalculationFile puts $num_step

    # TIME PER LOADSTEP
    GiD_WriteCalculationFile puts {$$TIME_PER_LOADSTEP}
    set time_step_path {string(container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="time_per_calculation_step"]/@v)}
    set time_step [$stageNode selectNodes $time_step_path]
    GiD_WriteCalculationFile puts $time_step

    # TOTAL TIME
    GiD_WriteCalculationFile puts {$$TOTAL_TIME}
    GiD_WriteCalculationFile puts $total_time
    
    # COURANT NUMBER
    GiD_WriteCalculationFile puts {$$COURANT_NUMBER}
    set courant_path {string(container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="Courant_number"]/@v)}
    set courant [$stageNode selectNodes $courant_path]
    GiD_WriteCalculationFile puts $courant

    # COMPUTATION METHOD
    GiD_WriteCalculationFile puts {$$COMPUTATION_METHOD}
    set method_path {string(container[@n="Calculation_Data"]/value[@n="COMPUTATION_METHOD"]/@v)}
    set method [$stageNode selectNodes $method_path]
    if {$method == "MPM - material point integration"} {
      GiD_WriteCalculationFile puts "MPM-MP"
    } elseif {$method == "MPM - mixed integration"} {
      GiD_WriteCalculationFile puts "MPM-MIXED"
    } elseif {$method == "standard FEM"} {
      GiD_WriteCalculationFile puts "FEM"
    } elseif {$method == "Updated Lagrangian FEM"} {
      GiD_WriteCalculationFile puts "UL-FEM"}   

    # GRAVITY ACCELERATION
    GiD_WriteCalculationFile puts {$$GRAVITY_ACCELERATION}
    set gacceleration_path {string(container[@n="General_data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_acceleration"]/@v)}
    set gacceleration [$root selectNodes $gacceleration_path]
    GiD_WriteCalculationFile puts $gacceleration
        
    # GRAVITY VECTOR
    GiD_WriteCalculationFile puts {$$GRAVITY_VECTOR}
    set gvector_x_path {string(container[@n="General_data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_x"]/@v)}
    set gvector_x [$root selectNodes $gvector_x_path]
    set gvector_y_path {string(container[@n="General_data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_y"]/@v)}
    set gvector_y [$root selectNodes $gvector_y_path]
    set gvector_z_path {string(container[@n="General_data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_z"]/@v)}
    set gvector_z [$root selectNodes $gvector_z_path]
    GiD_WriteCalculationFile puts [= "%s %s %s" $gvector_x $gvector_y $gvector_z]

    # GRAVITY LOAD
    GiD_WriteCalculationFile puts {$$GRAVITY_LOAD}
    set gload_path {string(container[@n="Calculation_Data"]/value[@n="GRAVITY_LOAD"]/@v)}
    set gload [$stageNode selectNodes $gload_path]
    set gmult_init_path {string(container[@n="Calculation_Data"]/value[@n="gravity_multiplier_initial"]/@v)}
    set gmult_init [$stageNode selectNodes $gmult_init_path]
    set gmult_fin_path {string(container[@n="Calculation_Data"]/value[@n="gravity_multiplier_final"]/@v)}
    set gmult_fin [$stageNode selectNodes $gmult_fin_path]
    if {$gload == "do not apply gravity load"} {
      set gload_type "off"    
    } elseif {$gload == "apply gravity load - linear"} {
      set gload_type "linear"             
    } elseif {$gload == "apply gravity load - stepwise"} {
      set gload_type "step"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $gload_type $gmult_init $gmult_fin]

    # SOLID TRACTION_A
    GiD_WriteCalculationFile puts {$$SOLID_TRACTION}
    set SolidTraction_A_path {string(container[@n="Calculation_Data"]/value[@n="SOLID_TRACTION_A"]/@v)}
    set SolidTraction_A [$stageNode selectNodes $SolidTraction_A_path]
    set SolidTraction_A_init_path {string(container[@n="Calculation_Data"]/value[@n="solid_traction_A_multiplier_initial"]/@v)}
    set SolidTraction_A_init [$stageNode selectNodes $SolidTraction_A_init_path]
    set SolidTraction_A_fin_path {string(container[@n="Calculation_Data"]/value[@n="solid_traction_A_multiplier_final"]/@v)}
    set SolidTraction_A_fin [$stageNode selectNodes $SolidTraction_A_fin_path]
    if {$SolidTraction_A == "do not apply solid traction"} {
      set SolidTraction_A_type "off"
    } elseif {$SolidTraction_A == "apply solid traction - linear"} {
      set SolidTraction_A_type "linear"       
    } elseif {$SolidTraction_A == "apply solid traction - stepwise"} {
      set SolidTraction_A_type "step"}  
    GiD_WriteCalculationFile puts [= "%s %s %s" $SolidTraction_A_type $SolidTraction_A_init $SolidTraction_A_fin]

    # SOLID TRACTION B
    GiD_WriteCalculationFile puts {$$SOLID_TRACTION_B}
    set SolidTraction_B_path {string(container[@n="Calculation_Data"]/value[@n="SOLID_TRACTION_B"]/@v)}
    set SolidTraction_B [$stageNode selectNodes $SolidTraction_B_path]
    set SolidTraction_B_init_path {string(container[@n="Calculation_Data"]/value[@n="solid_traction_B_multiplier_initial"]/@v)}
    set SolidTraction_B_init [$stageNode selectNodes $SolidTraction_B_init_path]
    set SolidTraction_B_fin_path {string(container[@n="Calculation_Data"]/value[@n="solid_traction_B_multiplier_final"]/@v)}
    set SolidTraction_B_fin [$stageNode selectNodes $SolidTraction_B_fin_path]
    if {$SolidTraction_B == "do not apply solid traction"} {
      set SolidTraction_B_type "off"
    } elseif {$SolidTraction_B == "apply solid traction - linear"} {
      set SolidTraction_B_type "linear" 
    } elseif {$SolidTraction_B == "apply solid traction - stepwise"} {
      set SolidTraction_B_type "step"}  
    GiD_WriteCalculationFile puts [= "%s %s %s" $SolidTraction_B_type $SolidTraction_B_init $SolidTraction_B_fin]
    
    # LIQUID PRESSURE_A
    GiD_WriteCalculationFile puts {$$LIQUID_PRESSURE}
    set LiquidPressure_A_path {string(container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_A"]/@v)}
    set LiquidPressure_A [$stageNode selectNodes $LiquidPressure_A_path]
    set LiquidPressure_A_init_path {string(container[@n="Calculation_Data"]/value[@n="liquid_pressure_A_multiplier_initial"]/@v)}
    set LiquidPressure_A_init [$stageNode selectNodes $LiquidPressure_A_init_path]
    set LiquidPressure_A_fin_path {string(container[@n="Calculation_Data"]/value[@n="liquid_pressure_A_multiplier_initial"]/@v)}
    set LiquidPressure_A_fin [$stageNode selectNodes $LiquidPressure_A_fin_path]
    if {$LiquidPressure_A == "do not apply liquid pressure"} {
      set LiquidPressure_A_type "off"
    } elseif {$LiquidPressure_A == "apply liquid pressure - linear"} {
      set LiquidPressure_A_type "linear"        
    } elseif {$LiquidPressure_A == "apply liquid pressure - stepwise"} {
      set LiquidPressure_A_type "step"}
      GiD_WriteCalculationFile puts [= "%s %s %s" $LiquidPressure_A_type $LiquidPressure_A_init $LiquidPressure_A_fin]  
    
    # LIQUID PRESSURE B
    GiD_WriteCalculationFile puts {$$LIQUID_PRESSURE_B}
    set LiquidPressure_B_path {string(container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_B"]/@v)}
    set LiquidPressure_B [$stageNode selectNodes $LiquidPressure_B_path]
    set LiquidPressure_B_init_path {string(container[@n="Calculation_Data"]/value[@n="liquid_pressure_B_multiplier_initial"]/@v)}
    set LiquidPressure_B_init [$stageNode selectNodes $LiquidPressure_B_init_path]
    set LiquidPressure_B_fin_path {string(container[@n="Calculation_Data"]/value[@n="liquid_pressure_B_multiplier_initial"]/@v)}
    set LiquidPressure_B_fin [$stageNode selectNodes $LiquidPressure_B_fin_path]
    if {$LiquidPressure_B == "do not apply liquid pressure"} {
      set LiquidPressure_B_type "off"
    } elseif {$LiquidPressure_B == "apply liquid pressure - linear"} {
      set LiquidPressure_B_type "linear"        
    } elseif {$LiquidPressure_B == "apply liquid pressure - stepwise"} {
      set LiquidPressure_B_type "step"}
      GiD_WriteCalculationFile puts [= "%s %s %s" $LiquidPressure_B_type $LiquidPressure_B_init $LiquidPressure_B_fin]

    # PRESCRIBED VELOCITY
    GiD_WriteCalculationFile puts {$$PRESCRIBED_VELOCITY}
    set PrescribedVelocity_path {string(container[@n="Calculation_Data"]/value[@n="PRESCRIBED_VELOCITY"]/@v)}
    set PrescribedVelocity [$stageNode selectNodes $PrescribedVelocity_path]
    set PrescribedVelocity_init_path {string(container[@n="Calculation_Data"]/value[@n="velocity_multiplier_initial"]/@v)}
    set PrescribedVelocity_init [$stageNode selectNodes $PrescribedVelocity_init_path]
    set PrescribedVelocity_fin_path {string(container[@n="Calculation_Data"]/value[@n="velocity_multiplier_final"]/@v)}
    set PrescribedVelocity_fin [$stageNode selectNodes $PrescribedVelocity_fin_path]
    if {$PrescribedVelocity == "do not apply prescribed velocity"} {
      set PrescribedVelocity_type "off"
    } elseif {$PrescribedVelocity == "apply prescribed velocity - linear"} {
      set PrescribedVelocity_type "linear"        
    } elseif {$PrescribedVelocity == "apply prescribed velocity - stepwise"} {
      set PrescribedVelocity_type "step"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $PrescribedVelocity_type $PrescribedVelocity_init $PrescribedVelocity_fin]
    
    # HYDRAULIC HEAD 
    GiD_WriteCalculationFile puts {$$HYDRAULIC_HEAD}
    set xp_hydraulic_head {container[@n="Calculation_Data"]/value[@n="HYDRAULIC_HEAD"]}
    set vNode [$stageNode selectNodes $xp_hydraulic_head]
    set v [get_domnode_attribute $vNode v]
    if { $v } {
        set HydraulicHead_type "file"
    } else {
        set HydraulicHead_type "off"
    }
    set HydraulicHead_value "0.0"
    GiD_WriteCalculationFile puts [= "%s %s" $HydraulicHead_type $HydraulicHead_value]

    # APPLY SEEPAGE FACE
    GiD_WriteCalculationFile puts {$$APPLY_SEEPAGE_FACE}
    set xp_seepage_face {container[@n="Calculation_Data"]/value[@n="SEEPAGE_FACE"]}
    set vNode [$stageNode selectNodes $xp_seepage_face]
    set v [get_domnode_attribute $vNode v]
    set SeepageFace_ID $v
    GiD_WriteCalculationFile puts $SeepageFace_ID
    
    # APPLY INFILTRATION
    GiD_WriteCalculationFile puts {$$APPLY_INFILTRATION}
    set xp_infiltration {container[@n="Calculation_Data"]/value[@n="INFILTRATION"]}
    set vNode [$stageNode selectNodes $xp_infiltration]
    set v [get_domnode_attribute $vNode v]
    set Infiltration_ID $v
    GiD_WriteCalculationFile puts $Infiltration_ID   

    # QUASI-STATIC CONVERGENCE
    GiD_WriteCalculationFile puts {$$QUASISTATIC_CONVERGENCE}
    set QuasiStatic_path {string(container[@n="Calculation_Data"]/value[@n="QUASI-STATIC_CONVERGENCE"]/@v)}
    set QuasiStatic [$stageNode selectNodes $QuasiStatic_path]
    if {$QuasiStatic == "do not apply convergence criteria"} {
        set QuasiStatic_ID "0"
    } elseif {$QuasiStatic == "apply convergence criteria"} {
        set QuasiStatic_ID "1"}
    GiD_WriteCalculationFile puts $QuasiStatic_ID    

    # TOLERATED ERROR SOLID, LIQUID, MAXIMUM TIME STEP 
    if {$QuasiStatic_ID == "1"} {      
      set SolidEnergy_path {string(container[@n="Calculation_Data"]/value[@n="tolerated_error_solid_energy"]/@v)}
      set SolidEnergy [$stageNode selectNodes $SolidEnergy_path]     
      set SolidForce_path {string(container[@n="Calculation_Data"]/value[@n="tolerated_error_solid_force"]/@v)}
      set SolidForce [$stageNode selectNodes $SolidForce_path]
      GiD_WriteCalculationFile puts {$$TOLERATED_ERROR_SOLID}      
      GiD_WriteCalculationFile puts [= "%s %s" $SolidEnergy $SolidForce]
      
      set LiquidEnergy_path {string(container[@n="Calculation_Data"]/value[@n="tolerated_error_liquid_energy"]/@v)}  
      set LiquidEnergy [$stageNode selectNodes $LiquidEnergy_path]  
      set LiquidForce_path {string(container[@n="Calculation_Data"]/value[@n="tolerated_error_liquid_force"]/@v)}  
      set LiquidForce [$stageNode selectNodes $LiquidForce_path]
      GiD_WriteCalculationFile puts {$$TOLERATED_ERROR_LIQUID}
      GiD_WriteCalculationFile puts [= "%s %s" $LiquidEnergy $LiquidForce]  
      
      set MaxTimeStep_path {string(container[@n="Calculation_Data"]/value[@n="maximum_time_steps"]/@v)}
      set MaxTimeStep [$stageNode selectNodes $MaxTimeStep_path]
      GiD_WriteCalculationFile puts {$$MAXIMUM_TIME_STEPS}
      GiD_WriteCalculationFile puts $MaxTimeStep
    }   

    # MASS SCALING
    GiD_WriteCalculationFile puts {$$MASS_SCALING}
    set MScaling_path {string(container[@n="Calculation_Data"]/value[@n="MASS_SCALING"]/@v)}
    set MScaling [$stageNode selectNodes $MScaling_path]
    set MScaling_value_path {string(container[@n="Calculation_Data"]/value[@n="mass_scaling_factor"]/@v)}
    set MScaling_value [$stageNode selectNodes $MScaling_value_path]
    if {$MScaling == "do not apply mass scaling"} {
      set MScaling_ID "0"
    } elseif {$MScaling == "apply mass scaling"} {
      set MScaling_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $MScaling_ID $MScaling_value]

    # HOMOGENEOUS LOCAL DAMPING
    GiD_WriteCalculationFile puts {$$HOMOGENEOUS_LOCAL_DAMPING}
    set LocalDamp_path {string(container[@n="Calculation_Data"]/value[@n="HOMOGENEOUS_LOCAL_DAMPING"]/@v)}
    set LocalDamp [$stageNode selectNodes $LocalDamp_path]
    set LocalDamp_value_path {string(container[@n="Calculation_Data"]/value[@n="local_damping_coefficient"]/@v)}
    set LocalDamp_value [$stageNode selectNodes $LocalDamp_value_path]
    if {$LocalDamp == "do not apply homogeneous local damping"} {
      set LocalDamp_ID "0"
    } elseif {$LocalDamp == "apply homogeneous local damping"} {
      set LocalDamp_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $LocalDamp_ID $LocalDamp_value]

    # BULK VISCOSITY DAMPING
    GiD_WriteCalculationFile puts {$$BULK_VISCOSITY_DAMPING}
    set ViscDamp_path {string(container[@n="Calculation_Data"]/value[@n="BULK_VISCOSITY_DAMPING"]/@v)}
    set ViscDamp [$stageNode selectNodes $ViscDamp_path]
    set Linear_path {string(container[@n="Calculation_Data"]/value[@n="linear_coefficient"]/@v)}
    set Linear [$stageNode selectNodes $Linear_path]
    set Quadratic_path {string(container[@n="Calculation_Data"]/value[@n="quadratic_coefficient"]/@v)}
    set Quadratic [$stageNode selectNodes $Quadratic_path]   
    if {$ViscDamp == "do not apply viscosity damping"} {
      set ViscDamp_ID "0"
    } elseif {$ViscDamp == "apply viscosity damping"} {
      set ViscDamp_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $ViscDamp_ID $Linear $Quadratic]   
    
    # STRAIN SMOOTHING
    GiD_WriteCalculationFile puts {$$STRAIN_SMOOTHING}
    set StrainSmoothing_path {string(container[@n="Calculation_Data"]/value[@n="STRAIN_SMOOTHING"]/@v)}
    set StrainSmoothing [$stageNode selectNodes $StrainSmoothing_path]
    if {$StrainSmoothing == "do not apply strain smoothing"} {
      set StrainSmoothing_ID "0"
    } elseif {$StrainSmoothing == "apply strain smoothing"} {
      set StrainSmoothing_ID "1"}
    GiD_WriteCalculationFile puts $StrainSmoothing_ID
    
    # SMOOTHENING LIQUID PRESSURE INCREMENT
    GiD_WriteCalculationFile puts {$$APPLY_SMOOTHENING_LIQUID_PRESSURE_INCREMENT}
    set LiquidPressureSmoothing_path {string(container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_INCREMENT_SMOOTHING"]/@v)}
    set LiquidPressureSmoothing [$stageNode selectNodes $LiquidPressureSmoothing_path]
    if {$LiquidPressureSmoothing == "do not apply liquid pressure increment smoothing"} {
      set LiquidPressureSmoothing_ID "0"
    } elseif {$LiquidPressureSmoothing == "apply liquid pressure increment smoothing"} {
      set LiquidPressureSmoothing_ID "1"}
    GiD_WriteCalculationFile puts $LiquidPressureSmoothing_ID
    
    # FIX SOLID SKELETON
    GiD_WriteCalculationFile puts {$$FIX_SOLID_SKELETON}
    set FixSolid_path {string(container[@n="Calculation_Data"]/value[@n="FIX_SOLID_SKELETON"]/@v)}
    set FixSolid [$stageNode selectNodes $FixSolid_path]
    if {$FixSolid == "Solid Skeleton Free"} {
      set FixSolid_ID "0"
    } elseif {$FixSolid == "Solid Skeleton Fixed"} {
      set FixSolid_ID "1"}
    GiD_WriteCalculationFile puts $FixSolid_ID

    # CONTACT FORMULATION
    GiD_WriteCalculationFile puts {$$CONTACT_FORMULATION}
    set Contact_pathNode [$stageNode selectNodes {container[@n="Calculation_Data"]/value[@n="CONTACT"]}]
    set Contact_ID [get_domnode_attribute $Contact_pathNode v]
    GiD_WriteCalculationFile puts $Contact_ID
    if {$Contact_ID == "1"} {
        set NormalCorrection_path {string(container[@n="Contact"]/value[@n="node_normal_correction"]/@v)}
        set NormalCorrection [$stageNode selectNodes $NormalCorrection_path]
        GiD_WriteCalculationFile puts {$$APPLY_CONTACT_NORMAL_CORRECTION}
        if {$NormalCorrection == "do not apply"} {
            set NormalCorrection_ID "0"
        } else {set NormalCorrection_ID "1"}
        GiD_WriteCalculationFile puts $NormalCorrection_ID
        GiD_WriteCalculationFile puts {$$NUMBER_OF_CORRECTED_NORMALS}
        set NumNormal_path {string(container[@n="Contact"]/value[@n="number_of_normals"]/@v)}
        set NumNormal [$stageNode selectNodes $NumNormal_path]
        if {$NumNormal == "none"} {
            GiD_WriteCalculationFile puts "0"
        } else {GiD_WriteCalculationFile puts $NumNormal}
        GiD_WriteCalculationFile puts {$$NODE_NORMAL_DATA}
        
        # Normal vector $icount       
        if {$NumNormal == "none"} {
            GiD_WriteCalculationFile puts "0"
        } 
        if {$NumNormal >= "1"} {
            Anura3D::WriteNormalVector $stageNode 1
        } 
        if {$NumNormal >= "2"} {
            Anura3D::WriteNormalVector $stageNode 2           
        } 
        if {$NumNormal >= "3"} {
            Anura3D::WriteNormalVector $stageNode 3          
        } 
        if {$NumNormal >= "4"} {
            Anura3D::WriteNormalVector $stageNode 4            
        } 
        if {$NumNormal >= "5"} {
            Anura3D::WriteNormalVector $stageNode 5            
        } 
        if {$NumNormal >= "6"} {
            Anura3D::WriteNormalVector $stageNode 6            
        } 
        if {$NumNormal >= "7"} {
            Anura3D::WriteNormalVector $stageNode 7            
        } 
        if {$NumNormal >= "8"} {
            Anura3D::WriteNormalVector $stageNode 8            
        } 
        if {$NumNormal >= "9"} {
            Anura3D::WriteNormalVector $stageNode 9            
        } 
        if {$NumNormal >= "10"} {
            Anura3D::WriteNormalVector $stageNode 10            
        }
    }
        
    # INITIAL WATER PRESSURE
    GiD_WriteCalculationFile puts {$$INITIAL_WATER_PRESSURE}
    set InitialWP_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/value[@n="INITIAL_WATER_PRESSURE"]/@v)}
    set InitialWP [$stageNode selectNodes $InitialWP_path]
    GiD_WriteCalculationFile puts $InitialWP
        
        # MATERIAL UPDATE
        GiD_WriteCalculationFile puts {$$APPLY_MATERIAL_UPDATE}
        set MaterialUpdate_path {string(value[@n="update_materials"]/@v)}
    set MaterialUpdate [$stageNode selectNodes $MaterialUpdate_path]
    if {$MaterialUpdate == "do not apply"} {
      set MaterialUpdate_ID "0"
    } elseif {$MaterialUpdate == "apply"} {
      set MaterialUpdate_ID "1"}
    GiD_WriteCalculationFile puts $MaterialUpdate_ID
    
    # POROSITY UPDATE
    GiD_WriteCalculationFile puts {$$APPLY_POROSITY_UPDATE}
    set PorosityUpdate_path {string(container[@n="Calculation_Data"]/value[@n="POROSITY_UPDATE"]/@v)}
    set PorosityUpdate [$stageNode selectNodes $PorosityUpdate_path]
    if {$PorosityUpdate == "do not apply porosity update"} {
      set PorosityUpdate_ID "0"
    } elseif {$PorosityUpdate == "apply porosity update"} {
      set PorosityUpdate_ID "1"}
    GiD_WriteCalculationFile puts $PorosityUpdate_ID
    
    # INITIAL VELOCITY
    GiD_WriteCalculationFile puts {$$INITIAL_VELOCITY}
    set xpInitialVelocity {container[@n="Calculation_Data"]/value[@n="INITIAL_VELOCITY"]}
    set InitialVelocity_IDNode [$stageNode selectNodes $xpInitialVelocity]
    set InitialVelocity_ID [get_domnode_attribute $InitialVelocity_IDNode v]                
    GiD_WriteCalculationFile puts $InitialVelocity_ID
    
    # RESET DISPLACEMENTS
    GiD_WriteCalculationFile puts {$$RESET_DISPLACEMENTS}
    set ResetDispl_path {string(container[@n="Calculation_Data"]/value[@n="RESET_DISPLACEMENTS"]/@v)}
    set ResetDispl [$stageNode selectNodes $ResetDispl_path]
    if {$ResetDispl == "do not reset displacements"} {
      set ResetDispl_ID "0"
    } elseif {$ResetDispl == "reset displacements"} {
      set ResetDispl_ID "1"}
    GiD_WriteCalculationFile puts $ResetDispl_ID
    
    # K0 PROCEDURE
    GiD_WriteCalculationFile puts {$$K0_PROCEDURE}
    set k0_procedure_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/value[@n="Apply_K0"]/@v)}
    set k0_procedure [$stageNode selectNodes $k0_procedure_path]
    set k0_procedure_ID "0"
    if {$k0_procedure == "apply"} {set k0_procedure_ID "1"}
    GiD_WriteCalculationFile puts $k0_procedure_ID
    if {$k0_procedure_ID == "1"} {      
        # Initial vertical load
        GiD_WriteCalculationFile puts {$$INITIAL_VERTICAL_LOAD_K0}
        set InitVerticalLoad_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/value[@n="initial_vertical_load"]/@v)}
        set InitVerticalLoad [$stageNode selectNodes $InitVerticalLoad_path]
        GiD_WriteCalculationFile puts $InitVerticalLoad
        # Max suction    
        set MaxSuction_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Not-horizontal"]/value[@n="MAX_SUCTION_AT_SOIL_SURFACE"]/@v)}
        set MaxSuction [$stageNode selectNodes $MaxSuction_path]
        set Suction_value_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Not-horizontal"]/value[@n="suction"]/@v)}
        set Suction_value [$stageNode selectNodes $Suction_value_path]
        if {$MaxSuction == "specify"} {
            GiD_WriteCalculationFile puts {$$K0_MAX_SUCTION}        
            GiD_WriteCalculationFile puts $Suction_value}
        # Surface elevation         
        set surf_elevation_path {string(//container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="soil_surface"]/@v)}
        set surf_elevation [$stageNode selectNodes $surf_elevation_path]        
        GiD_WriteCalculationFile puts {$$SURFACE_ELEVATION}
        GiD_WriteCalculationFile puts $surf_elevation
        # Number of layers
        set NumLayers_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="NUMBER_OF_LAYERS"]/@v)}
        set NumLayers [$stageNode selectNodes $NumLayers_path]
        GiD_WriteCalculationFile puts {$$NUMBER_SOIL_LAYERS}
        if {$NumLayers == "none"} {  
            GiD_WriteCalculationFile puts "0"
        } else {GiD_WriteCalculationFile puts $NumLayers}
        # Thickness of soil layers
        if {($NumLayers >= "1") && ($NumLayers != "none")} {
            GiD_WriteCalculationFile puts {$$THICKNESS_SOIL_LAYERS}
            set Layer1_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_1"]/@v)}
            set Layer1 [$stageNode selectNodes $Layer1_path]
            set Layer2_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_2"]/@v)}
            set Layer2 [$stageNode selectNodes $Layer2_path]
            set Layer3_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_3"]/@v)}
            set Layer3 [$stageNode selectNodes $Layer3_path]
            set Layer4_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_4"]/@v)}
            set Layer4 [$stageNode selectNodes $Layer4_path]
            set Layer5_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_5"]/@v)}
            set Layer5 [$stageNode selectNodes $Layer5_path]
            set Layer6_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_6"]/@v)}
            set Layer6 [$stageNode selectNodes $Layer6_path]
            set Layer7_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_7"]/@v)}
            set Layer7 [$stageNode selectNodes $Layer7_path]
            set Layer8_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_8"]/@v)}
            set Layer8 [$stageNode selectNodes $Layer8_path]
            set Layer9_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_9"]/@v)}
            set Layer9 [$stageNode selectNodes $Layer9_path]
            set Layer10_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/container[@n="Horizontal"]/value[@n="layer_thickness_10"]/@v)}
            set Layer10 [$stageNode selectNodes $Layer10_path]
            if {$NumLayers == "1"} {
                GiD_WriteCalculationFile puts $Layer1
            } elseif {$NumLayers == "2"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
            } elseif {$NumLayers == "3"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
            } elseif {$NumLayers == "4"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
            } elseif {$NumLayers == "5"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
            } elseif {$NumLayers == "6"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
                GiD_WriteCalculationFile puts $Layer6           
            } elseif {$NumLayers == "7"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
                GiD_WriteCalculationFile puts $Layer6
                GiD_WriteCalculationFile puts $Layer7        
            } elseif {$NumLayers == "8"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
                GiD_WriteCalculationFile puts $Layer6
                GiD_WriteCalculationFile puts $Layer7
                GiD_WriteCalculationFile puts $Layer8        
            } elseif {$NumLayers == "9"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
                GiD_WriteCalculationFile puts $Layer6
                GiD_WriteCalculationFile puts $Layer7
                GiD_WriteCalculationFile puts $Layer8
                GiD_WriteCalculationFile puts $Layer9 
            } elseif {$NumLayers == "10"} {
                GiD_WriteCalculationFile puts $Layer1
                GiD_WriteCalculationFile puts $Layer2
                GiD_WriteCalculationFile puts $Layer3
                GiD_WriteCalculationFile puts $Layer4
                GiD_WriteCalculationFile puts $Layer5
                GiD_WriteCalculationFile puts $Layer6
                GiD_WriteCalculationFile puts $Layer7
                GiD_WriteCalculationFile puts $Layer8
                GiD_WriteCalculationFile puts $Layer9
                GiD_WriteCalculationFile puts $Layer10 }  
        }    
    }
    
    # STRESS INITIALIZATION FROM FILE
    GiD_WriteCalculationFile puts {$$MP_STRESS_INITIALIZATION_FROMFILE}
    set StressFromFile_path {string(container[@n="Initial_cond"]/container[@n="Stress_initialization"]/value[@n="FROM_FILE"]/@v)}
    set StressFromFile [$stageNode selectNodes $StressFromFile_path]
    set StressFromFile_ID "0"
    if {$StressFromFile == "assign stresses from external file"} {set StressFromFile_ID "1"}
    GiD_WriteCalculationFile puts $StressFromFile_ID
        
    # EXCAVATION
    set num_tot 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {condition[@n="Solid_Excavation"]/group[@ov=%s]} $ov_type]
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
        set num_tot [expr $gEntities_num + $num_tot] 
    }
    if {$num_tot != 0} {
        GiD_WriteCalculationFile puts {$$EXCAVATION}    
        GiD_WriteCalculationFile puts $num_tot 
    }
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]  
        set gEntities_id [GiD_EntitiesGroups get $gName $ov_type]
        for {set i 0} {$i < $gEntities_num } {incr i} {
            set id_entity [lindex $gEntities_id $i]
            set FirstStep [$gNode selectNodes {string(value[@n="First_step"]/@v)}]
            set LastStep [$gNode selectNodes {string(value[@n="Last_step"]/@v)}]
            GiD_WriteCalculationFile puts [= "%s %s %s" $id_entity $FirstStep $LastStep]} 
    }

    # CONSTRUCTION
    set num_tot 0   
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
        set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set ov_type "volume"
    }
    set xp [format_xpath {container[@n="Construction"]/condition[@n="Solid_Construction"]/group[@ov=%s]} $ov_type]
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
        set num_tot [expr $gEntities_num + $num_tot] 
    }
    if {$num_tot != 0} {
        GiD_WriteCalculationFile puts {$$CONSTRUCTION}    
        GiD_WriteCalculationFile puts $num_tot 
    }
    foreach gNode [$stageNode selectNodes $xp] {
        set gName [get_domnode_attribute $gNode n]
        set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]  
        set gEntities_id [GiD_EntitiesGroups get $gName $ov_type]
        for {set i 0} {$i < $gEntities_num } {incr i} {
            set id_entity [lindex $gEntities_id $i]
            set FirstStep [$gNode selectNodes {string(value[@n="First_step"]/@v)}]
            set LastStep [$gNode selectNodes {string(value[@n="Last_step"]/@v)}]
            GiD_WriteCalculationFile puts [= "%s %s %s" $id_entity $FirstStep $LastStep]} 
    }

    # FILL EMPTY ELEMENTS
    if {$num_tot != 0} {
    GiD_WriteCalculationFile puts {$$FILL_EMPTY_ELEMENTS}
    set FillEmptyElem_path {string(container[@n="Construction"]/value[@n="Fill_empty_elements"]/@v)}
    set FillEmptyElem [$stageNode selectNodes $FillEmptyElem_path]
    if {$FillEmptyElem == "do not fill empty elements"} {
      set FillEmptyElem_ID "0"
    } elseif {$FillEmptyElem == "fill empty elements"} {
      set FillEmptyElem_ID "1"}
    GiD_WriteCalculationFile puts $FillEmptyElem_ID
    }

   # SUBMERGED CALCULATION
    GiD_WriteCalculationFile puts {$$SUBMERGED_CALCULATION}
    set Submerged_path {string(container[@n="Calculation_Data"]/value[@n="SUBMERGED_CALCULATION"]/@v)}
    set Submerged [$stageNode selectNodes $Submerged_path]
    set Submerged_steps_path {string(container[@n="Calculation_Data"]/value[@n="number_of_initialisation_steps"]/@v)}
    set Submerged_steps [$stageNode selectNodes $Submerged_steps_path]
    set Submerged_ID "0"
    if {$Submerged == "apply submerged calculation"} {set Submerged_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $Submerged_ID $Submerged_steps]
    
    # OBJECTIVE STRESS
    GiD_WriteCalculationFile puts {$$APPLY_OBJECTIVE_STRESS}
    set ObjectiveStress_path {string(container[@n="Calculation_Data"]/value[@n="OBJECTIVE_STRESS"]/@v)}
    set ObjectiveStress [$stageNode selectNodes $ObjectiveStress_path]
    set ObjectiveStress_ID "0"
    if {$ObjectiveStress == "apply objective stress"} {set ObjectiveStress_ID "1"}
    GiD_WriteCalculationFile puts $ObjectiveStress_ID
    
    # DEGREE OF FILLING 
    GiD_WriteCalculationFile puts {$$DEGREE_OF_FILLING}
    set DegreeOfFilling_path {string(container[@n="Calculation_Data"]/value[@n="degree_of_filling"]/@v)}
    set DegreeOfFilling [$stageNode selectNodes $DegreeOfFilling_path]
    GiD_WriteCalculationFile puts $DegreeOfFilling  
    
    # NUMBER OF ACTIVE ELEMENTS           
    GiD_WriteCalculationFile puts {$$NUMBER_OF_ACTIVE_ELEMENTS}
    GiD_WriteCalculationFile puts $number_of_active_elements
    
    # DOUBLE-POINT FORMULATION
    set DoublePoint_path {string(container[@n="Calculation_Data"]/value[@n="DOUBLE-POINT_FORMULATION"]/@v)}
    set DoublePoint [$stageNode selectNodes $DoublePoint_path]
    if {$DoublePoint == "define additional parameters"} {
      # Max porosity
      set MaxPorosity_path {string(container[@n="Calculation_Data"]/value[@n="maximum_porosity"]/@v)}
      set MaxPorosity [$stageNode selectNodes $MaxPorosity_path]
      GiD_WriteCalculationFile puts {$$MAXIMUM_POROSITY}
      GiD_WriteCalculationFile puts $MaxPorosity
      # Ergun constants
      set ErgunConst_path {string(container[@n="Calculation_Data"]/value[@n="Ergun_constants"]/@v)}
      set ErgunConst [$stageNode selectNodes $ErgunConst_path]
      GiD_WriteCalculationFile puts {$$ERGUN_CONSTANTS}
      GiD_WriteCalculationFile puts $ErgunConst 
      # Permeability update
      set Permeability_path {string(container[@n="Calculation_Data"]/value[@n="permeability_definition"]/@v)}
      set Permeability [$stageNode selectNodes $Permeability_path]
      GiD_WriteCalculationFile puts {$$PERMEABILITY_UPDATE}
          if {$Permeability == "please choose"} {
            set answer [tk_messageBox -title "Anura3D - Generate Anura3D Files" -message "INPUT ERROR: For the double-point formulation the permeability update has to be chosen." -icon warning -type okcancel]
      } elseif {$Permeability == "constant permeability"} {     
            GiD_WriteCalculationFile puts "constant_permeability"
            GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY}
            set IntrinsicPermeability_path {string(container[@n="Calculation_Data"]/value[@n="intrinsic_permeability"]/@v)}
            set IntrinsicPermeability [$stageNode selectNodes $IntrinsicPermeability_path]
            GiD_WriteCalculationFile puts $IntrinsicPermeability    
      } elseif {$Permeability == "update permeability Darcy"} {
            GiD_WriteCalculationFile puts "Darcy_update"  
            GiD_WriteCalculationFile puts {$$GRAIN_SIZE_DIAMETER}
            set GrainSizeDiameter_path {string(container[@n="Calculation_Data"]/value[@n="grain_size_diameter"]/@v)}
            set GrainSizeDiameter [$stageNode selectNodes $GrainSizeDiameter_path]
            GiD_WriteCalculationFile puts $GrainSizeDiameter  
      } elseif {$Permeability == "update permeability Ergun"} {
            GiD_WriteCalculationFile puts "Ergun_update"
            GiD_WriteCalculationFile puts {$$GRAIN_SIZE_DIAMETER}
            set GrainSizeDiameter_path {string(container[@n="Calculation_Data"]/value[@n="grain_size_diameter"]/@v)}
            set GrainSizeDiameter [$stageNode selectNodes $GrainSizeDiameter_path]
            GiD_WriteCalculationFile puts $GrainSizeDiameter
      }   
      # Elevation liquid surface k0
      GiD_WriteCalculationFile puts {$$LIQUID_SURFACE}
      set LiquidSurface_path {string(container[@n="Calculation_Data"]/value[@n="elevation_liquid_surface_for_K0"]/@v)}
      set LiquidSurface [$stageNode selectNodes $LiquidSurface_path]
      GiD_WriteCalculationFile puts $LiquidSurface
      # Strain smoothing liquid 
      GiD_WriteCalculationFile puts {$$APPLY_STRAIN_SMOOTHING_LIQUID_TWOLAYERFORM}
      set StrainSmoothingLiquid_path {string(container[@n="Calculation_Data"]/value[@n="strain_smoothing_liquid"]/@v)}
      set StrainSmoothingLiquid [$stageNode selectNodes $StrainSmoothingLiquid_path]
      GiD_WriteCalculationFile puts $StrainSmoothingLiquid        
      # No tensile stress liquid
      GiD_WriteCalculationFile puts {$$NO_TENSILE_STRESS_LIQUID_MP_WITH_LIQUID_STATUS}
      set NoTensileStress_path {string(container[@n="Calculation_Data"]/value[@n="no_tensile_stress_liquid"]/@v)}
      set NoTensileStress [$stageNode selectNodes $NoTensileStress_path]
      GiD_WriteCalculationFile puts $NoTensileStress       
      # Detect free surface
      GiD_WriteCalculationFile puts {$$DETECT_FREE_SURFACE}
      set DetectFreeSurface_path {string(container[@n="Calculation_Data"]/value[@n="detect_free_surface"]/@v)}
      set DetectFreeSurface [$stageNode selectNodes $DetectFreeSurface_path]
      GiD_WriteCalculationFile puts $DetectFreeSurface
    }    
    
    # OUTPUT VISUALIZATION DATA    
    GiD_WriteCalculationFile puts {$$VISUALIZATION_OPTION}
    set paraview [$root selectNodes {string(container[@n="General_data"]/container[@n="POSTPROCESS_VISUALIZATION_SOFTWARE"]/value[@n="Paraview"]/@v)}]
    set GiD_ASCII [$root selectNodes {string(container[@n="General_data"]/container[@n="POSTPROCESS_VISUALIZATION_SOFTWARE"]/value[@n="GiD_ASCII"]/@v)}]
    set GiD_Binary [$root selectNodes {string(container[@n="General_data"]/container[@n="POSTPROCESS_VISUALIZATION_SOFTWARE"]/value[@n="GiD_Binary"]/@v)}]
   
    if { $paraview && !$GiD_ASCII && !$GiD_Binary } {
        set output "Paraview"       
    } 
    if { $GiD_ASCII } {
        if { !$paraview } {
            set output "GiD-ASCII"
        } else {
            set output "Paraview-GiD"
        }         
    } 
    if { $GiD_Binary } {
        if { !$paraview } {
            set output "GiD-Binary"
        } else {
            set output "Paraview-GiD"
        }                 
    }
    GiD_WriteCalculationFile puts $output 

    # OUTPUT MP DATA
    if { $icount_stage != "0"} { 
        GiD_WriteCalculationFile puts {$$OUTPUT_NUMBER_OF_MATERIAL_POINTS}
        GiD_WriteCalculationFile puts $out_num_material_points
		if { $out_num_material_points > "0"} {   
			GiD_WriteCalculationFile puts {$$OUTPUT_MATERIAL_POINTS}
			set i 0
			for {set i 0} {$i < $out_num_material_points } {incr i} {
				set mat_points_id [lindex $out_material_points  $i]
				GiD_WriteCalculationFile puts $mat_points_id} 
		}
	}
   
    GiD_WriteCalculationFile puts {$$END}
    GiD_WriteCalculationFile end
}

proc Anura3D::WriteNormalVector { stageNode i } {    
    set format "" 
    set value "normal_vector_$i"
    set xp_normvect [format_xpath {container[@n="Contact"]/condition[@n=%s]/group} $value]                    
    foreach gNode [$stageNode selectNodes $xp_normvect] {
        set n [$gNode @n]
        dict set format $n "%d\n"
    }
    if { [dict size $format] } {
        set err [catch {
                GiD_WriteCalculationFile nodes -unique -return $format} numnodes]       
        if { $err } {          
            snit_messageBox -parent .gid \
                -message [= "Error when writing nodes for normal_vector $i (%s)" $numnodes]             
            return "isnotok"                        
        } else {                                       
            GiD_WriteCalculationFile nodes -unique $format
        }
    }
}
    