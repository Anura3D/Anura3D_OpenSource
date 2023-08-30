#      CPS file
#      print data in the .dat calculation file (instead of a classic .bas template)


proc Anura3D::WriteCalculationFile_CPS { filename } {
  variable current_xml_root
  GiD_WriteCalculationFile init $filename
    set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
    set exe_name [GiD_Info Project ProblemType]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data
    customlib::SetBaseRoot $root
    set current_xml_root $root

  GiD_WriteCalculationFile puts "### Anura3D_2023 ###"

    # NUMBER OF LOADSTEP
    GiD_WriteCalculationFile puts {$$NUMBER_OF_LOADSTEPS}
    set num_step_path {string(//container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="number_of_calculation_steps"]/@v)}
    set num_step [$root selectNodes $num_step_path]
    GiD_WriteCalculationFile puts $num_step

    # TIME PER LOADSTEP
    GiD_WriteCalculationFile puts {$$TIME_PER_LOADSTEP}
    set time_step_path {string(//container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="time_per_calculation_step"]/@v)}
    set time_step [$root selectNodes $time_step_path]
    GiD_WriteCalculationFile puts $time_step

    # TOTAL TIME
    GiD_WriteCalculationFile puts {$$TOTAL_TIME}
    GiD_WriteCalculationFile puts "0.0"
    
    # COURANT NUMBER
    GiD_WriteCalculationFile puts {$$COURANT_NUMBER}
    set courant_path {string(//container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/value[@n="Courant_number"]/@v)}
    set courant [$root selectNodes $courant_path]
    GiD_WriteCalculationFile puts $courant

    # COMPUTATION METHOD
    GiD_WriteCalculationFile puts {$$COMPUTATION_METHOD}
    set method_path {string(//container[@n="Calculation_Data"]/value[@n="COMPUTATION_METHOD"]/@v)}
    set method [$root selectNodes $method_path]
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
    set gacceleration_path {string(//container[@n="Calculation_Data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_acceleration"]/@v)}
    set gacceleration [$root selectNodes $gacceleration_path]
    GiD_WriteCalculationFile puts $gacceleration
	
    # GRAVITY VECTOR
    GiD_WriteCalculationFile puts {$$GRAVITY_VECTOR}
    set gvector_x_path {string(//container[@n="Calculation_Data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_x"]/@v)}
    set gvector_x [$root selectNodes $gvector_x_path]
    set gvector_y_path {string(//container[@n="Calculation_Data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_y"]/@v)}
    set gvector_y [$root selectNodes $gvector_y_path]
    set gvector_z_path {string(//container[@n="Calculation_Data"]/container[@n="GRAVITY_DATA"]/value[@n="gravity_vector_z"]/@v)}
    set gvector_z [$root selectNodes $gvector_z_path]
    GiD_WriteCalculationFile puts [= "%s %s %s" $gvector_x $gvector_y $gvector_z]

    # GRAVITY LOAD
    GiD_WriteCalculationFile puts {$$GRAVITY_LOAD}
    set gload_path {string(//container[@n="Calculation_Data"]/value[@n="GRAVITY_LOAD"]/@v)}
    set gload [$root selectNodes $gload_path]
    set gmult_init_path {string(//container[@n="Calculation_Data"]/value[@n="gravity_multiplier_initial"]/@v)}
    set gmult_init [$root selectNodes $gmult_init_path]
    set gmult_fin_path {string(//container[@n="Calculation_Data"]/value[@n="gravity_multiplier_final"]/@v)}
    set gmult_fin [$root selectNodes $gmult_fin_path]
    if {$gload == "do not apply gravity load"} {
      set gload_type "off"    
    } elseif {$gload == "apply gravity load - linear"} {
      set gload_type "linear"             
    } elseif {$gload == "apply gravity load - stepwise"} {
      set gload_type "step"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $gload_type $gmult_init $gmult_fin]

    # SOLID TRACTION_A
    GiD_WriteCalculationFile puts {$$SOLID_TRACTION}
    set SolidTraction_A_path {string(//container[@n="Calculation_Data"]/value[@n="SOLID_TRACTION_A"]/@v)}
    set SolidTraction_A [$root selectNodes $SolidTraction_A_path]
    set SolidTraction_A_init_path {string(//container[@n="Calculation_Data"]/value[@n="solid_traction_A_multiplier_initial"]/@v)}
    set SolidTraction_A_init [$root selectNodes $SolidTraction_A_init_path]
    set SolidTraction_A_fin_path {string(//container[@n="Calculation_Data"]/value[@n="solid_traction_A_multiplier_final"]/@v)}
    set SolidTraction_A_fin [$root selectNodes $SolidTraction_A_fin_path]
    if {$SolidTraction_A == "do not apply solid traction"} {
      set SolidTraction_A_type "off"
    } elseif {$SolidTraction_A == "apply solid traction - linear"} {
      set SolidTraction_A_type "linear"       
    } elseif {$SolidTraction_A == "apply solid traction - stepwise"} {
      set SolidTraction_A_type "step"}  
    GiD_WriteCalculationFile puts [= "%s %s %s" $SolidTraction_A_type $SolidTraction_A_init $SolidTraction_A_fin]

    # SOLID TRACTION B
    GiD_WriteCalculationFile puts {$$SOLID_TRACTION_B}
    set SolidTraction_B_path {string(//container[@n="Calculation_Data"]/value[@n="SOLID_TRACTION_B"]/@v)}
    set SolidTraction_B [$root selectNodes $SolidTraction_B_path]
    set SolidTraction_B_init_path {string(//container[@n="Calculation_Data"]/value[@n="solid_traction_B_multiplier_initial"]/@v)}
    set SolidTraction_B_init [$root selectNodes $SolidTraction_B_init_path]
    set SolidTraction_B_fin_path {string(//container[@n="Calculation_Data"]/value[@n="solid_traction_B_multiplier_final"]/@v)}
    set SolidTraction_B_fin [$root selectNodes $SolidTraction_B_fin_path]
    if {$SolidTraction_B == "do not apply solid traction"} {
      set SolidTraction_B_type "off"
    } elseif {$SolidTraction_B == "apply solid traction - linear"} {
      set SolidTraction_B_type "linear" 
    } elseif {$SolidTraction_B == "apply solid traction - stepwise"} {
      set SolidTraction_B_type "step"}  
    GiD_WriteCalculationFile puts [= "%s %s %s" $SolidTraction_B_type $SolidTraction_B_init $SolidTraction_B_fin]
    
    # LIQUID PRESSURE_A
    GiD_WriteCalculationFile puts {$$LIQUID_PRESSURE}
    set LiquidPressure_A_path {string(//container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_A"]/@v)}
    set LiquidPressure_A [$root selectNodes $LiquidPressure_A_path]
    set LiquidPressure_A_init_path {string(//container[@n="Calculation_Data"]/value[@n="liquid_pressure_A_multiplier_initial"]/@v)}
    set LiquidPressure_A_init [$root selectNodes $LiquidPressure_A_init_path]
    set LiquidPressure_A_fin_path {string(//container[@n="Calculation_Data"]/value[@n="liquid_pressure_A_multiplier_initial"]/@v)}
    set LiquidPressure_A_fin [$root selectNodes $LiquidPressure_A_fin_path]
    if {$LiquidPressure_A == "do not apply liquid pressure"} {
      set LiquidPressure_A_type "off"
    } elseif {$LiquidPressure_A == "apply liquid pressure - linear"} {
      set LiquidPressure_A_type "linear"        
    } elseif {$LiquidPressure_A == "apply liquid pressure - stepwise"} {
      set LiquidPressure_A_type "step"}
      GiD_WriteCalculationFile puts [= "%s %s %s" $LiquidPressure_A_type $LiquidPressure_A_init $LiquidPressure_A_fin]  
    
    # LIQUID PRESSURE B
    GiD_WriteCalculationFile puts {$$LIQUID_PRESSURE_B}
    set LiquidPressure_B_path {string(//container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_B"]/@v)}
    set LiquidPressure_B [$root selectNodes $LiquidPressure_B_path]
    set LiquidPressure_B_init_path {string(//container[@n="Calculation_Data"]/value[@n="liquid_pressure_B_multiplier_initial"]/@v)}
    set LiquidPressure_B_init [$root selectNodes $LiquidPressure_B_init_path]
    set LiquidPressure_B_fin_path {string(//container[@n="Calculation_Data"]/value[@n="liquid_pressure_B_multiplier_initial"]/@v)}
    set LiquidPressure_B_fin [$root selectNodes $LiquidPressure_B_fin_path]
    if {$LiquidPressure_B == "do not apply liquid pressure"} {
      set LiquidPressure_B_type "off"
    } elseif {$LiquidPressure_B == "apply liquid pressure - linear"} {
      set LiquidPressure_B_type "linear"        
    } elseif {$LiquidPressure_B == "apply liquid pressure - stepwise"} {
      set LiquidPressure_B_type "step"}
      GiD_WriteCalculationFile puts [= "%s %s %s" $LiquidPressure_B_type $LiquidPressure_B_init $LiquidPressure_B_fin]

    # PRESCRIBED VELOCITY
    GiD_WriteCalculationFile puts {$$PRESCRIBED_VELOCITY}
    set PrescribedVelocity_path {string(//container[@n="Calculation_Data"]/value[@n="PRESCRIBED_VELOCITY"]/@v)}
    set PrescribedVelocity [$root selectNodes $PrescribedVelocity_path]
    set PrescribedVelocity_init_path {string(//container[@n="Calculation_Data"]/value[@n="velocity_multiplier_initial"]/@v)}
    set PrescribedVelocity_init [$root selectNodes $PrescribedVelocity_init_path]
    set PrescribedVelocity_fin_path {string(//container[@n="Calculation_Data"]/value[@n="velocity_multiplier_final"]/@v)}
    set PrescribedVelocity_fin [$root selectNodes $PrescribedVelocity_fin_path]
    if {$PrescribedVelocity == "do not apply prescribed velocity"} {
      set PrescribedVelocity_type "off"
    } elseif {$PrescribedVelocity == "apply prescribed velocity - linear"} {
      set PrescribedVelocity_type "linear"        
    } elseif {$PrescribedVelocity == "apply prescribed velocity - stepwise"} {
      set PrescribedVelocity_type "step"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $PrescribedVelocity_type $PrescribedVelocity_init $PrescribedVelocity_fin]

    # HYDRAULIC HEAD
    GiD_WriteCalculationFile puts {$$HYDRAULIC_HEAD}
    set HydraulicHead_path {string(//container[@n="Calculation_Data"]/value[@n="HYDRAULIC_HEAD"]/@v)}
    set HydraulicHead [$root selectNodes $HydraulicHead_path]
    set HydraulicHead_value "0.0"
    if {$HydraulicHead == "do not apply hydraulic head"} {
      set HydraulicHead_type "off"
    } elseif {$HydraulicHead == "apply hydraulic head"} {
      set HydraulicHead_type "file"}    
    GiD_WriteCalculationFile puts [= "%s %s" $HydraulicHead_type $HydraulicHead_value]
    

    # APPLY SEEPAGE FACE
    GiD_WriteCalculationFile puts {$$APPLY_SEEPAGE_FACE}
    set SeepageFace_path {string(//container[@n="Calculation_Data"]/value[@n="SEEPAGE_FACE"]/@v)}
    set SeepageFace [$root selectNodes $SeepageFace_path]
    if {$SeepageFace == "do not apply seepage face"} {
	set SeepageFace_ID "0"
    } elseif {$SeepageFace == "apply seepage face"} {
	set SeepageFace_ID "1"}
    GiD_WriteCalculationFile puts $SeepageFace_ID        

    # APPLY INFILTRATION
    GiD_WriteCalculationFile puts {$$APPLY_INFILTRATION}
    set Infiltration_path {string(//container[@n="Calculation_Data"]/value[@n="INFILTRATION"]/@v)}
    set Infiltration [$root selectNodes $Infiltration_path]
    if {$Infiltration == "do not apply infiltration"} {
	set Infiltration_ID "0"
    } elseif {$Infiltration == "apply infiltration"} {
	set Infiltration_ID "1"}
    GiD_WriteCalculationFile puts $Infiltration_ID   

    # QUASI-STATIC CONVERGENCE
    GiD_WriteCalculationFile puts {$$QUASISTATIC_CONVERGENCE}
    set QuasiStatic_path {string(//container[@n="Calculation_Data"]/value[@n="QUASI-STATIC_CONVERGENCE"]/@v)}
    set QuasiStatic [$root selectNodes $QuasiStatic_path]
    if {$QuasiStatic == "do not apply convergence criteria"} {
	set QuasiStatic_ID "0"
    } elseif {$QuasiStatic == "apply convergence criteria"} {
	set QuasiStatic_ID "1"}
    GiD_WriteCalculationFile puts $QuasiStatic_ID    

    # TOLERATED ERROR SOLID, LIQUID, MAXIMUM TIME STEP 
    if {$QuasiStatic_ID == "1"} {      
      set SolidEnergy_path {string(//container[@n="Calculation_Data"]/value[@n="tolerated_error_solid_energy"]/@v)}
      set SolidEnergy [$root selectNodes $SolidEnergy_path]     
      set SolidForce_path {string(//container[@n="Calculation_Data"]/value[@n="tolerated_error_solid_force"]/@v)}
      set SolidForce [$root selectNodes $SolidForce_path]
      GiD_WriteCalculationFile puts {$$TOLERATED_ERROR_SOLID}      
      GiD_WriteCalculationFile puts [= "%s %s" $SolidEnergy $SolidForce]
      
      set LiquidEnergy_path {string(//container[@n="Calculation_Data"]/value[@n="tolerated_error_liquid_energy"]/@v)}  
      set LiquidEnergy [$root selectNodes $LiquidEnergy_path]  
      set LiquidForce_path {string(//container[@n="Calculation_Data"]/value[@n="tolerated_error_liquid_force"]/@v)}  
      set LiquidForce [$root selectNodes $LiquidForce_path]
      GiD_WriteCalculationFile puts {$$TOLERATED_ERROR_LIQUID}
      GiD_WriteCalculationFile puts [= "%s %s" $LiquidEnergy $LiquidForce]  
      
      set MaxTimeStep_path {string(//container[@n="Calculation_Data"]/value[@n="maximum_time_steps"]/@v)}
      set MaxTimeStep [$root selectNodes $MaxTimeStep_path]
      GiD_WriteCalculationFile puts {$$MAXIMUM_TIME_STEPS}
      GiD_WriteCalculationFile puts $MaxTimeStep
    }   

    # MASS SCALING
    GiD_WriteCalculationFile puts {$$MASS_SCALING}
    set MScaling_path {string(//container[@n="Calculation_Data"]/value[@n="MASS_SCALING"]/@v)}
    set MScaling [$root selectNodes $MScaling_path]
    set MScaling_value_path {string(//container[@n="Calculation_Data"]/value[@n="mass_scaling_factor"]/@v)}
    set MScaling_value [$root selectNodes $MScaling_value_path]
    if {$MScaling == "do not apply mass scaling"} {
      set MScaling_ID "0"
    } elseif {$MScaling == "apply mass scaling"} {
      set MScaling_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $MScaling_ID $MScaling_value]

    # HOMOGENEOUS LOCAL DAMPING
    GiD_WriteCalculationFile puts {$$HOMOGENEOUS_LOCAL_DAMPING}
    set LocalDamp_path {string(//container[@n="Calculation_Data"]/value[@n="HOMOGENEOUS_LOCAL_DAMPING"]/@v)}
    set LocalDamp [$root selectNodes $LocalDamp_path]
    set LocalDamp_value_path {string(//container[@n="Calculation_Data"]/value[@n="local_damping_coefficient"]/@v)}
    set LocalDamp_value [$root selectNodes $LocalDamp_value_path]
    if {$LocalDamp == "do not apply homogeneous local damping"} {
      set LocalDamp_ID "0"
    } elseif {$LocalDamp == "apply homogeneous local damping"} {
      set LocalDamp_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $LocalDamp_ID $LocalDamp_value]

    # BULK VISCOSITY DAMPING
    GiD_WriteCalculationFile puts {$$BULK_VISCOSITY_DAMPING}
    set ViscDamp_path {string(//container[@n="Calculation_Data"]/value[@n="BULK_VISCOSITY_DAMPING"]/@v)}
    set ViscDamp [$root selectNodes $ViscDamp_path]
    set Linear_path {string(//container[@n="Calculation_Data"]/value[@n="linear_coefficient"]/@v)}
    set Linear [$root selectNodes $Linear_path]
    set Quadratic_path {string(//container[@n="Calculation_Data"]/value[@n="quadratic_coefficient"]/@v)}
    set Quadratic [$root selectNodes $Quadratic_path]   
    if {$ViscDamp == "do not apply viscosity damping"} {
      set ViscDamp_ID "0"
    } elseif {$ViscDamp == "apply viscosity damping"} {
      set ViscDamp_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $ViscDamp_ID $Linear $Quadratic]   
    
    # STRAIN SMOOTHING
    GiD_WriteCalculationFile puts {$$STRAIN_SMOOTHING}
    set StrainSmoothing_path {string(//container[@n="Calculation_Data"]/value[@n="STRAIN_SMOOTHING"]/@v)}
    set StrainSmoothing [$root selectNodes $StrainSmoothing_path]
    if {$StrainSmoothing == "do not apply strain smoothing"} {
      set StrainSmoothing_ID "0"
    } elseif {$StrainSmoothing == "apply strain smoothing"} {
      set StrainSmoothing_ID "1"}
    GiD_WriteCalculationFile puts $StrainSmoothing_ID
    
    # SMOOTHENING LIQUID PRESSURE INCREMENT
    GiD_WriteCalculationFile puts {$$APPLY_SMOOTHENING_LIQUID_PRESSURE_INCREMENT}
    set LiquidPressureSmoothing_path {string(//container[@n="Calculation_Data"]/value[@n="LIQUID_PRESSURE_INCREMENT_SMOOTHING"]/@v)}
    set LiquidPressureSmoothing [$root selectNodes $LiquidPressureSmoothing_path]
    if {$LiquidPressureSmoothing == "do not apply liquid pressure increment smoothing"} {
      set LiquidPressureSmoothing_ID "0"
    } elseif {$LiquidPressureSmoothing == "apply liquid pressure increment smoothing"} {
      set LiquidPressureSmoothing_ID "1"}
    GiD_WriteCalculationFile puts $LiquidPressureSmoothing_ID
    
    # FIX SOLID SKELETON
    GiD_WriteCalculationFile puts {$$FIX_SOLID_SKELETON}
    set FixSolid_path {string(//container[@n="Calculation_Data"]/value[@n="FIX_SOLID_SKELETON"]/@v)}
    set FixSolid [$root selectNodes $FixSolid_path]
    if {$FixSolid == "Solid Skeleton Free"} {
      set FixSolid_ID "0"
    } elseif {$FixSolid == "Solid Skeleton Fixed"} {
      set FixSolid_ID "1"}
    GiD_WriteCalculationFile puts $FixSolid_ID

    # CONTACT FORMULATION
    GiD_WriteCalculationFile puts {$$CONTACT_FORMULATION}
    set Contact_path {string(//container[@n="Calculation_Data"]/value[@n="CONTACT"]/@v)}
    set Contact [$root selectNodes $Contact_path]
    if {$Contact == "do not apply contact"} {
      set Contact_ID "0"
    } elseif {$Contact == "apply contact"} {
      set Contact_ID "1"}
	GiD_WriteCalculationFile puts $Contact_ID
	if {$Contact_ID == "1"} {
	set NormalCorrection_path {string(//container[@n="Calculation_Data"]/value[@n="node_normal_correction"]/@v)}
    set NormalCorrection [$root selectNodes $NormalCorrection_path]
	GiD_WriteCalculationFile puts {$$APPLY_CONTACT_NORMAL_CORRECTION}
	if {$NormalCorrection == "do not apply node normal correction"} {
	  set NormalCorrection_ID "0"
	} else {set NormalCorrection_ID "1"}
	GiD_WriteCalculationFile puts $NormalCorrection_ID
	GiD_WriteCalculationFile puts {$$NUMBER_OF_CORRECTED_NORMALS}
	set NumNormal_path {string(//container[@n="Calculation_Data"]/value[@n="number_of_normals"]/@v)}
    set NumNormal [$root selectNodes $NumNormal_path]
    if {$NumNormal == "none"} {
      GiD_WriteCalculationFile puts "0"
    } else {GiD_WriteCalculationFile puts $NumNormal}
    GiD_WriteCalculationFile puts {$$NODE_NORMAL_DATA}
    set Vector1_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_1"]/@v)}
    set Vector1 [$root selectNodes $Vector1_path]
	set Vector2_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_2"]/@v)}
    set Vector2 [$root selectNodes $Vector2_path]
	set Vector3_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_3"]/@v)}
    set Vector3 [$root selectNodes $Vector3_path]
	set Vector4_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_4"]/@v)}
    set Vector4 [$root selectNodes $Vector4_path]
	set Vector5_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_5"]/@v)}
    set Vector5 [$root selectNodes $Vector5_path]
	set Vector6_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_6"]/@v)}
    set Vector6 [$root selectNodes $Vector6_path]
	set Vector7_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_7"]/@v)}
    set Vector7 [$root selectNodes $Vector7_path]
	set Vector8_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_8"]/@v)}
    set Vector8 [$root selectNodes $Vector8_path]
	set Vector9_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_9"]/@v)}
    set Vector9 [$root selectNodes $Vector9_path]
	set Vector10_path {string(//container[@n="Calculation_Data"]/value[@n="normal_vector_10"]/@v)}
    set Vector10 [$root selectNodes $Vector10_path]	
	if {$NumNormal == "none"} {
      GiD_WriteCalculationFile puts "0"
    } elseif {$NumNormal == "1"} {
      GiD_WriteCalculationFile puts $Vector1
	} elseif {$NumNormal == "2"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	} elseif {$NumNormal == "3"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	} elseif {$NumNormal == "4"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	} elseif {$NumNormal == "5"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	} elseif {$NumNormal == "6"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	  GiD_WriteCalculationFile puts $Vector6	
	} elseif {$NumNormal == "7"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	  GiD_WriteCalculationFile puts $Vector6
	  GiD_WriteCalculationFile puts $Vector7
	} elseif {$NumNormal == "8"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	  GiD_WriteCalculationFile puts $Vector6
	  GiD_WriteCalculationFile puts $Vector7
	  GiD_WriteCalculationFile puts $Vector8	
	} elseif {$NumNormal == "9"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	  GiD_WriteCalculationFile puts $Vector6
	  GiD_WriteCalculationFile puts $Vector7
	  GiD_WriteCalculationFile puts $Vector8
	  GiD_WriteCalculationFile puts $Vector9	
	} elseif {$NumNormal == "10"} {
	  GiD_WriteCalculationFile puts $Vector1
	  GiD_WriteCalculationFile puts $Vector2
	  GiD_WriteCalculationFile puts $Vector3
	  GiD_WriteCalculationFile puts $Vector4
	  GiD_WriteCalculationFile puts $Vector5
	  GiD_WriteCalculationFile puts $Vector6
	  GiD_WriteCalculationFile puts $Vector7
	  GiD_WriteCalculationFile puts $Vector8
	  GiD_WriteCalculationFile puts $Vector9
	  GiD_WriteCalculationFile puts $Vector10 }
    }
	
    # INITIAL WATER PRESSURE
    GiD_WriteCalculationFile puts {$$INITIAL_WATER_PRESSURE}
    set InitialWP_path {string(//container[@n="Calculation_Data"]/value[@n="INITIAL_WATER_PRESSURE"]/@v)}
    set InitialWP [$root selectNodes $InitialWP_path]
    set InitialWP_value_path {string(//container[@n="Calculation_Data"]/value[@n="water_pressure"]/@v)}
    set InitialWP_value [$root selectNodes $InitialWP_value_path]
    GiD_WriteCalculationFile puts $InitialWP_value
    
    # POROSITY UPDATE
    GiD_WriteCalculationFile puts {$$APPLY_POROSITY_UPDATE}
    set PorosityUpdate_path {string(//container[@n="Calculation_Data"]/value[@n="POROSITY_UPDATE"]/@v)}
    set PorosityUpdate [$root selectNodes $PorosityUpdate_path]
    if {$PorosityUpdate == "do not apply porosity update"} {
      set PorosityUpdate_ID "0"
    } elseif {$PorosityUpdate == "apply porosity update"} {
      set PorosityUpdate_ID "1"}
    GiD_WriteCalculationFile puts $PorosityUpdate_ID
    
    # INITIAL VELOCITY
    GiD_WriteCalculationFile puts {$$INITIAL_VELOCITY}
    set InitialVelocity_path {string(//container[@n="Calculation_Data"]/value[@n="INITIAL_VELOCITY"]/@v)}
    set InitialVelocity [$root selectNodes $InitialVelocity_path]
    if {$InitialVelocity == "do not apply initial velocity"} {
      set InitialVelocity_ID "0"
    } elseif {$InitialVelocity == "apply initial velocity"} {
      set InitialVelocity_ID "1"}
    GiD_WriteCalculationFile puts $InitialVelocity_ID
    
    # RESET DISPLACEMENTS
    GiD_WriteCalculationFile puts {$$RESET_DISPLACEMENTS}
    set ResetDispl_path {string(//container[@n="Calculation_Data"]/value[@n="RESET_DISPLACEMENTS"]/@v)}
    set ResetDispl [$root selectNodes $ResetDispl_path]
    if {$ResetDispl == "do not reset displacements"} {
      set ResetDispl_ID "0"
    } elseif {$ResetDispl == "reset displacements"} {
      set ResetDispl_ID "1"}
    GiD_WriteCalculationFile puts $ResetDispl_ID
    
    # REMOVE FIXITIES
    GiD_WriteCalculationFile puts {$$REMOVE_FIXITIES}
    set RemoveFixitySolid_path {string(//container[@n="Calculation_Data"]/value[@n="REMOVE_FIXITIES"]/value[@n="solid"]/@v)}
    set RemoveFixitySolid [$root selectNodes $RemoveFixitySolid_path]
    set RemoveFixitySolid_ID "0"
    if {$RemoveFixitySolid == "remove fixities"} {set RemoveFixitySolid_ID "1"}    
    set RemoveFixityLiquid_path {string(//container[@n="Calculation_Data"]/value[@n="REMOVE_FIXITIES"]/value[@n="liquid"]/@v)}
    set RemoveFixityLiquid [$root selectNodes $RemoveFixityLiquid_path]
    set RemoveFixityLiquid_ID "0"
    if {$RemoveFixityLiquid == "remove fixities"} {set RemoveFixityLiquid_ID "1"}   
    set RemoveFixityGas_path {string(//container[@n="Calculation_Data"]/value[@n="REMOVE_FIXITIES"]/value[@n="gas"]/@v)}
    set RemoveFixityGas [$root selectNodes $RemoveFixityGas_path]
    set RemoveFixityGas_ID "0"
    if {$RemoveFixityGas == "remove fixities"} {set RemoveFixityGas_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s %s" $RemoveFixitySolid_ID $RemoveFixityLiquid_ID $RemoveFixityGas_ID]
      
    # K0 PROCEDURE
    GiD_WriteCalculationFile puts {$$K0_PROCEDURE}
    set k0_procedure_path {string(//container[@n="Calculation_Data"]/value[@n="K0-PROCEDURE"]/@v)}
    set k0_procedure [$root selectNodes $k0_procedure_path]
    set k0_procedure_ID "0"
    if {$k0_procedure == "apply K0-procedure"} {set k0_procedure_ID "1"}
    GiD_WriteCalculationFile puts $k0_procedure_ID
    if {$k0_procedure_ID == "1"} {      
    # Surface elevation
    GiD_WriteCalculationFile puts {$$SURFACE_ELEVATION}
    set SurfaceElevation_path {string(//container[@n="Calculation_Data"]/value[@n="soil_surface"]/@v)}
    set SurfaceElevation [$root selectNodes $SurfaceElevation_path]
    GiD_WriteCalculationFile puts $SurfaceElevation
    # Initial vertical load
    GiD_WriteCalculationFile puts {$$INITIAL_VERTICAL_LOAD_K0}
    set InitVerticalLoad_path {string(//container[@n="Calculation_Data"]/value[@n="initial_vertical_load"]/@v)}
    set InitVerticalLoad [$root selectNodes $InitVerticalLoad_path]
    GiD_WriteCalculationFile puts $InitVerticalLoad
    # Max suction    
    set MaxSuction_path {string(//container[@n="Calculation_Data"]/value[@n="MAX_SUCTION_AT_SOIL_SURFACE"]/@v)}
    set MaxSuction [$root selectNodes $MaxSuction_path]
    set Suction_value_path {string(//container[@n="Calculation_Data"]/value[@n="suction"]/@v)}
    set Suction_value [$root selectNodes $Suction_value_path]
    if {$MaxSuction == "specify"} {
    GiD_WriteCalculationFile puts {$$K0_MAX_SUCTION}	
    GiD_WriteCalculationFile puts $Suction_value
	}
    # Number of layers
    set NumLayers_path {string(//container[@n="Calculation_Data"]/value[@n="NUMBER_OF_LAYERS"]/@v)}
    set NumLayers [$root selectNodes $NumLayers_path]
    GiD_WriteCalculationFile puts {$$NUMBER_SOIL_LAYERS}
	if {$NumLayers == "none"} {  
    GiD_WriteCalculationFile puts "0"
	} else {GiD_WriteCalculationFile puts $NumLayers}
    # Thickness of soil layers
	if {($NumLayers >= "1") && ($NumLayers != "none")} {
    GiD_WriteCalculationFile puts {$THICKNESS_SOIL_LAYERS}
    set Layer1_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_1"]/@v)}
    set Layer1 [$root selectNodes $Layer1_path]
    set Layer2_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_2"]/@v)}
    set Layer2 [$root selectNodes $Layer2_path]
    set Layer3_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_3"]/@v)}
    set Layer3 [$root selectNodes $Layer3_path]
    set Layer4_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_4"]/@v)}
    set Layer4 [$root selectNodes $Layer4_path]
    set Layer5_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_5"]/@v)}
    set Layer5 [$root selectNodes $Layer5_path]
    set Layer6_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_6"]/@v)}
    set Layer6 [$root selectNodes $Layer6_path]
    set Layer7_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_7"]/@v)}
    set Layer7 [$root selectNodes $Layer7_path]
    set Layer8_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_8"]/@v)}
    set Layer8 [$root selectNodes $Layer8_path]
    set Layer9_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_9"]/@v)}
    set Layer9 [$root selectNodes $Layer9_path]
    set Layer10_path {string(//container[@n="Calculation_Data"]/value[@n="layer_thickness_10"]/@v)}
    set Layer10 [$root selectNodes $Layer10_path]
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
    
    # EXCAVATION
    set num_tot 0
    set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
    set dim_type [$current_xml_root selectNodes $dim_path]
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	set ov_type "surface"
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	set ov_type "volume"
    }
    set xp [format_xpath {condition[@n="Solid_Excavation"]/group[@ov=%s]} $ov_type]
    foreach gNode [$root selectNodes $xp] {
	set gName [get_domnode_attribute $gNode n]
	set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
	set num_tot [expr $gEntities_num + $num_tot] 
    }
    if {$num_tot != 0} {
	GiD_WriteCalculationFile puts {$$EXCAVATION}    
	GiD_WriteCalculationFile puts $num_tot 
    }
    foreach gNode [$root selectNodes $xp] {
	set gName [get_domnode_attribute $gNode n]
	set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]  
	set gEntities_id [GiD_EntitiesGroups get $gName $ov_type]
	for {set i 0} {$i < $gEntities_num } {incr i} {
	    set id_entity [lindex $gEntities_id $i]
	    set FirstStep [$gNode selectNodes {string(value[@n="First_step"]/@v)}]
	    set LastStep [$gNode selectNodes {string(value[@n="Last_step"]/@v)}]
	    GiD_WriteCalculationFile puts [= "%s %s %s" $id_entity $FirstStep $LastStep]} 
    }

    # SUBMERGED CALCULATION
    GiD_WriteCalculationFile puts {$$SUBMERGED_CALCULATION}
    set Submerged_path {string(//container[@n="Calculation_Data"]/value[@n="SUBMERGED_CALCULATION"]/@v)}
    set Submerged [$root selectNodes $Submerged_path]
    set Submerged_steps_path {string(//container[@n="Calculation_Data"]/value[@n="number_of_initialisation_steps"]/@v)}
    set Submerged_steps [$root selectNodes $Submerged_steps_path]
    set Submerged_ID "0"
    if {$Submerged == "apply submerged calculation"} {set Submerged_ID "1"}
    GiD_WriteCalculationFile puts [= "%s %s" $Submerged_ID $Submerged_steps]
    
    # OBJECTIVE STRESS
    GiD_WriteCalculationFile puts {$$APPLY_OBJECTIVE_STRESS}
    set ObjectiveStress_path {string(//container[@n="Calculation_Data"]/value[@n="OBJECTIVE_STRESS"]/@v)}
    set ObjectiveStress [$root selectNodes $ObjectiveStress_path]
    set ObjectiveStress_ID "0"
    if {$ObjectiveStress == "apply objective stress"} {set ObjectiveStress_ID "1"}
    GiD_WriteCalculationFile puts $ObjectiveStress_ID
    
    # DEGREE OF FILLING 
    GiD_WriteCalculationFile puts {$$DEGREE_OF_FILLING}
    set DegreeOfFilling_path {string(//container[@n="Calculation_Data"]/value[@n="degree_of_filling"]/@v)}
    set DegreeOfFilling [$root selectNodes $DegreeOfFilling_path]
    GiD_WriteCalculationFile puts $DegreeOfFilling  
    
    # NUMBER OF ACTIVE ELEMENTS
	set i 0
	set layer_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
    set layer_type [$current_xml_root selectNodes $layer_path]
	if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {
	  set ov_type "surface"
	  
	  if {$layer_type == "Double_point"} {
	    set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Double-point"]/group} $ov_type]
	  } else {
	    set xp [format_xpath {container[@n="MPspecification"]/condition[@n="2D_Single-point"]/group} $ov_type]
	  }        
	} elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
	  set ov_type "volume"
	  if {$layer_type == "Double_point"} {
	    set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Double-point"]/group} $ov_type]
	  } else {
	  set xp [format_xpath {container[@n="MPspecification"]/condition[@n="3D_Single-point"]/group} $ov_type]
	  }        
	}        
	foreach gNode [$root selectNodes $xp] {       		
	  set l_material [$gNode selectNodes {string(value[@n="material"]/@v)}]
	  set list_group [$gNode @n]
	  set elements_id [GiD_EntitiesGroups get $list_group elements]
	  set num_elems [objarray length $elements_id]
	  set i [expr $i + $num_elems]
	}	
	GiD_WriteCalculationFile puts {$$NUMBER_OF_ACTIVE_ELEMENTS}
	GiD_WriteCalculationFile puts $i
    
    # DOUBLE-POINT FORMULATION
    set DoublePoint_path {string(//container[@n="Calculation_Data"]/value[@n="DOUBLE-POINT_FORMULATION"]/@v)}
    set DoublePoint [$root selectNodes $DoublePoint_path]
    if {$DoublePoint == "define additional parameters"} {
      # Max porosity
      set MaxPorosity_path {string(//container[@n="Calculation_Data"]/value[@n="maximum_porosity"]/@v)}
      set MaxPorosity [$root selectNodes $MaxPorosity_path]
      GiD_WriteCalculationFile puts {$$MAXIMUM_POROSITY}
      GiD_WriteCalculationFile puts $MaxPorosity
      # Ergun constants
      set ErgunConst_path {string(//container[@n="Calculation_Data"]/value[@n="Ergun_constants"]/@v)}
      set ErgunConst [$root selectNodes $ErgunConst_path]
      GiD_WriteCalculationFile puts {$$ERGUN_CONSTANTS}
      GiD_WriteCalculationFile puts $ErgunConst 
      # Permeability update
      set Permeability_path {string(//container[@n="Calculation_Data"]/value[@n="permeability_definition"]/@v)}
      set Permeability [$root selectNodes $Permeability_path]
      GiD_WriteCalculationFile puts {$$PERMEABILITY_UPDATE}
	  if {$Permeability == "please choose"} {
	    set answer [tk_messageBox -title "Anura3D - Generate Anura3D Files" -message "INPUT ERROR: For the double-point formulation the permeability update has to be chosen." -icon warning -type okcancel]
      } elseif {$Permeability == "constant permeability"} {     
	    GiD_WriteCalculationFile puts "constant_permeability"
	    GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY}
	    set IntrinsicPermeability_path {string(//container[@n="Calculation_Data"]/value[@n="intrinsic_permeability"]/@v)}
	    set IntrinsicPermeability [$root selectNodes $IntrinsicPermeability_path]
	    GiD_WriteCalculationFile puts $IntrinsicPermeability    
      } elseif {$Permeability == "update permeability Darcy"} {
	    GiD_WriteCalculationFile puts "Darcy_update"  
	    GiD_WriteCalculationFile puts {$$GRAIN_SIZE_DIAMETER}
	    set GrainSizeDiameter_path {string(//container[@n="Calculation_Data"]/value[@n="grain_size_diameter"]/@v)}
	    set GrainSizeDiameter [$root selectNodes $GrainSizeDiameter_path]
	    GiD_WriteCalculationFile puts $GrainSizeDiameter  
      } elseif {$Permeability == "update permeability Ergun"} {
	    GiD_WriteCalculationFile puts "Ergun_update"
	    GiD_WriteCalculationFile puts {$$GRAIN_SIZE_DIAMETER}
	    set GrainSizeDiameter_path {string(//container[@n="Calculation_Data"]/value[@n="grain_size_diameter"]/@v)}
	    set GrainSizeDiameter [$root selectNodes $GrainSizeDiameter_path]
	    GiD_WriteCalculationFile puts $GrainSizeDiameter
      }   
      # Elevation liquid surface k0
      GiD_WriteCalculationFile puts {$$LIQUID_SURFACE}
      set LiquidSurface_path {string(//container[@n="Calculation_Data"]/value[@n="elevation_liquid_surface_for_K0"]/@v)}
      set LiquidSurface [$root selectNodes $LiquidSurface_path]
      GiD_WriteCalculationFile puts $LiquidSurface
      # Strain smoothing liquid 
      GiD_WriteCalculationFile puts {$$APPLY_STRAIN_SMOOTHING_LIQUID_TWOLAYERFORM}
      set StrainSmoothingLiquid_path {string(//container[@n="Calculation_Data"]/value[@n="strain_smoothing_liquid"]/@v)}
      set StrainSmoothingLiquid [$root selectNodes $StrainSmoothingLiquid_path]
      GiD_WriteCalculationFile puts $StrainSmoothingLiquid        
      # No tensile stress liquid
      GiD_WriteCalculationFile puts {$$NO_TENSILE_STRESS_LIQUID_MP_WITH_LIQUID_STATUS}
      set NoTensileStress_path {string(//container[@n="Calculation_Data"]/value[@n="no_tensile_stress_liquid"]/@v)}
      set NoTensileStress [$root selectNodes $NoTensileStress_path]
      GiD_WriteCalculationFile puts $NoTensileStress       
      # Detect free surface
      GiD_WriteCalculationFile puts {$$DETECT_FREE_SURFACE}
      set DetectFreeSurface_path {string(//container[@n="Calculation_Data"]/value[@n="detect_free_surface"]/@v)}
      set DetectFreeSurface [$root selectNodes $DetectFreeSurface_path]
      GiD_WriteCalculationFile puts $DetectFreeSurface
    }    
    
    # OUTPUT VISUALIZATION DATA
    GiD_WriteCalculationFile puts {$$VISUALIZATION_OPTION}
    set OutputVisualization_path {string(//container[@n="Calculation_Data"]/value[@n="POSTPROCESS_VISUALIZATION_OPTIONS"]/@v)}
    set OutputVisualization [$root selectNodes $OutputVisualization_path]
    if {$OutputVisualization == "Paraview visualization"} {
	set output "Paraview"
    } elseif {$OutputVisualization == "GiD visualization ASCII format"} {
	set output "GiD-ASCII"
    } elseif {$OutputVisualization == "GiD visualization Binary format"} {
	set output "GiD-Binary"
    } elseif {$OutputVisualization == "GiD and Paraview"} {
	set output "Paraview-GiD"}   
    GiD_WriteCalculationFile puts $output
    
    # OUTPUT DATA
    GiD_WriteCalculationFile puts {$$OUTPUT_NUMBER_OF_MATERIAL_POINTS}
    set NumMP_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="number_of_material_points"]/@v)}
    set NumMP [$root selectNodes $NumMP_path]
    if {$NumMP == "none"} {set number_of_MP "0"} {set number_of_MP $NumMP}
    GiD_WriteCalculationFile puts $number_of_MP
    set MP1_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_1"]/@v)}
    set MP1 [$root selectNodes $MP1_path]
    set MP2_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_2"]/@v)}
    set MP2 [$root selectNodes $MP2_path]
    set MP3_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_3"]/@v)}
    set MP3 [$root selectNodes $MP3_path]
    set MP4_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_4"]/@v)}
    set MP4 [$root selectNodes $MP4_path]
    set MP5_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_5"]/@v)}
    set MP5 [$root selectNodes $MP5_path]
    set MP6_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_6"]/@v)}
    set MP6 [$root selectNodes $MP6_path]
    set MP7_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_7"]/@v)}
    set MP7 [$root selectNodes $MP7_path]
    set MP8_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_8"]/@v)}
    set MP8 [$root selectNodes $MP8_path]
    set MP9_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_9"]/@v)}
    set MP9 [$root selectNodes $MP9_path]
    set MP10_path {string(//container[@n="Calculation_Data"]/container[@n="OUTPUT_DATA"]/value[@n="material_point_10"]/@v)}
    set MP10 [$root selectNodes $MP10_path]
    GiD_WriteCalculationFile puts {$$OUTPUT_MATERIAL_POINTS}
    if {$NumMP == "none"} {
      GiD_WriteCalculationFile puts "0"
    } elseif {$NumMP == "1"} {
      GiD_WriteCalculationFile puts $MP1
    } elseif {$NumMP == "2"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2    
    } elseif {$NumMP == "3"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3    
    } elseif {$NumMP == "4"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4  
    } elseif {$NumMP == "5"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5   
    } elseif {$NumMP == "6"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5
      GiD_WriteCalculationFile puts $MP6  
    } elseif {$NumMP == "7"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5
      GiD_WriteCalculationFile puts $MP6
      GiD_WriteCalculationFile puts $MP7  
    } elseif {$NumMP == "8"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5
      GiD_WriteCalculationFile puts $MP6
      GiD_WriteCalculationFile puts $MP7
      GiD_WriteCalculationFile puts $MP8  
    } elseif {$NumMP == "9"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5
      GiD_WriteCalculationFile puts $MP6
      GiD_WriteCalculationFile puts $MP7
      GiD_WriteCalculationFile puts $MP8
      GiD_WriteCalculationFile puts $MP9  
    } elseif {$NumMP == "10"} {
      GiD_WriteCalculationFile puts $MP1
      GiD_WriteCalculationFile puts $MP2
      GiD_WriteCalculationFile puts $MP3
      GiD_WriteCalculationFile puts $MP4
      GiD_WriteCalculationFile puts $MP5
      GiD_WriteCalculationFile puts $MP6
      GiD_WriteCalculationFile puts $MP7
      GiD_WriteCalculationFile puts $MP8
      GiD_WriteCalculationFile puts $MP9
      GiD_WriteCalculationFile puts $MP10
    }    
  GiD_WriteCalculationFile puts {$$END}
GiD_WriteCalculationFile end
}
