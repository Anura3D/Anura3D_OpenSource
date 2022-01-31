### Anura3D_2021 ###
$$NUMBER_OF_LOADSTEPS
*GenData(number_of_calculation_steps_[-]__)
$$TIME_PER_LOADSTEP
*GenData(time_per_calculation_step_[s]__)
$$TOTAL_TIME
0.0
$$COURANT_NUMBER
*GenData(Courant_number_[-]__)
$$COMPUTATION_METHOD
*if((GenData(computation_ID,int)==100))
MPM-MIXED
*end if
*if((GenData(computation_ID,int)==101))
MPM-MP
*end if
*if((GenData(computation_ID,int)==102))
FEM
*end if
*if((GenData(computation_ID,int)==103))
UL-FEM
*end if
$$GRAVITY_ACCELERATION
*GenData(gravity_acceleration_[m/s2]__)
$$GRAVITY_VECTOR
*GenData(gravity_vector_x_[-]__) *GenData(gravity_vector_y_[-]__) *GenData(gravity_vector_z_[-]__)
$$GRAVITY_LOAD
*GenData(gravity_load_type)  *GenData(gravity_multiplier_initial_[-])  *GenData(gravity_multiplier_final_[-])
$$SOLID_TRACTION
*GenData(solid_traction_type)  *GenData(solid_traction_multiplier_initial_[-])  *GenData(solid_traction_multiplier_final_[-])
$$LIQUID_PRESSURE
*GenData(liquid_pressure_type)  *GenData(liquid_pressure_multiplier_initial_[-])  *GenData(liquid_pressure_multiplier_final_[-])
$$PRESCRIBED_VELOCITY
*GenData(prescribed_velocity_type)  *GenData(velocity_multiplier_initial_[-])  *GenData(velocity_multiplier_final_[-])
$$QUASISTATIC_CONVERGENCE
*GenData(quasi-static_convergence_ID)
$$TOLERATED_ERROR_SOLID
*GenData(tolerated_error_solid_energy_[-]__) *GenData(tolerated_error_solid_force_[-]__)
$$TOLERATED_ERROR_LIQUID
*GenData(tolerated_error_liquid_energy_[-]__) *GenData(tolerated_error_liquid_force_[-]__)
$$MAXIMUM_TIME_STEPS
*GenData(Maximum_time_steps_[-]__)
$$MASS_SCALING
*GenData(mass_scaling_ID) *GenData(mass_scaling_factor_[-]__)
$$HOMOGENEOUS_LOCAL_DAMPING
*GenData(damping_ID) *GenData(local_damping_coefficient_[-]__)
$$BULK_VISCOSITY_DAMPING
*GenData(bulk_viscosity_ID) *GenData(linear_coefficient_[-]__) *GenData(quadratic_coefficient_[-]__)
$$FIX_SOLID_SKELETON
*GenData(fix_solid_skeleton_ID)
$$STRAIN_SMOOTHING
*GenData(smoothing_ID)
$$CONTACT_FORMULATION
*GenData(contact_ID)
$$INITIAL_WATER_PRESSURE
*GenData(water_pressure_[kPa])
$$APPLY_POROSITY_UPDATE
*GenData(porosity_update_ID)
$$INITIAL_VELOCITY
*GenData(initialvelocity_ID)
$$RESET_DISPLACEMENTS
*GenData(reset_displacement_ID)
$$REMOVE_FIXITIES
*GenData(remove_fixities_solid_ID) *GenData(remove_fixities_liquid_ID) *GenData(remove_fixities_gas_ID)
*##### K0-procedure #####
*if((GenData(k0_procedure_ID,int)==1))
$$K0_PROCEDURE
*GenData(k0_procedure_ID)
$$SURFACE_ELEVATION
*GenData(soil_surface_[m]__)
$$INITIAL_VERTICAL_LOAD_K0
*GenData(initial_vertical_load_[kPa]__)
*end if
*##### excavation #####
$$EXCAVATION
*if(ndime==3)
*set cond 3D_-_Solid_Excavation *elems
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var tmp5=5
*set var tmp6=6
*set var tmp7=7
*set var tmp8=8
*set var tmp9=9
*set var tmp10=10
*set var nExcavation=0
*loop elems onlyincond
*set var VolumeNumber=cond(2,int)
*if(VolumeNumber==tmp1)
*set var tmp1=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp2)
*set var tmp2=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp3)
*set var tmp3=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp4)
*set var tmp4=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp5)
*set var tmp5=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp6)
*set var tmp6=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp7)
*set var tmp7=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp8)
*set var tmp8=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp9)
*set var tmp9=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp10)
*set var tmp10=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*end elems
*nExcavation
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var tmp5=5
*set var tmp6=6
*set var tmp7=7
*set var tmp8=8
*set var tmp9=9
*set var tmp10=10
*set var nExcavation=0
*loop elems onlyincond
*set var VolumeNumber=cond(2,int)
*if(VolumeNumber==tmp1)
*cond(2) *cond(3) *cond(4)
*set var tmp1=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp2)
*cond(2) *cond(3) *cond(4)
*set var tmp2=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp3)
*cond(2) *cond(3) *cond(4)
*set var tmp3=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp4)
*cond(2) *cond(3) *cond(4)
*set var tmp4=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp5)
*cond(2) *cond(3) *cond(4)
*set var tmp5=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp6)
*cond(2) *cond(3) *cond(4)
*set var tmp6=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp7)
*cond(2) *cond(3) *cond(4)
*set var tmp7=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp8)
*cond(2) *cond(3) *cond(4)
*set var tmp8=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp9)
*cond(2) *cond(3) *cond(4)
*set var tmp9=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp10)
*cond(2) *cond(3) *cond(4)
*set var tmp10=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*end elems
*end if
*if(ndime==2)
*set cond 2D_-_Solid_Excavation *elems
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var tmp5=5
*set var tmp6=6
*set var tmp7=7
*set var tmp8=8
*set var tmp9=9
*set var tmp10=10
*set var nExcavation=0
*loop elems onlyincond
*set var VolumeNumber=cond(2,int)
*if(VolumeNumber==tmp1)
*set var tmp1=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp2)
*set var tmp2=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp3)
*set var tmp3=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp4)
*set var tmp4=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp5)
*set var tmp5=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp6)
*set var tmp6=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp7)
*set var tmp7=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp8)
*set var tmp8=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp9)
*set var tmp9=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp10)
*set var tmp10=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*end elems
*nExcavation
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var tmp5=5
*set var tmp6=6
*set var tmp7=7
*set var tmp8=8
*set var tmp9=9
*set var tmp10=10
*set var nExcavation=0
*loop elems onlyincond
*set var VolumeNumber=cond(2,int)
*if(VolumeNumber==tmp1)
*cond(2) *cond(3) *cond(4)
*set var tmp1=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp2)
*cond(2) *cond(3) *cond(4)
*set var tmp2=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp3)
*cond(2) *cond(3) *cond(4)
*set var tmp3=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp4)
*cond(2) *cond(3) *cond(4)
*set var tmp4=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp5)
*cond(2) *cond(3) *cond(4)
*set var tmp5=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp6)
*cond(2) *cond(3) *cond(4)
*set var tmp6=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp7)
*cond(2) *cond(3) *cond(4)
*set var tmp7=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp8)
*cond(2) *cond(3) *cond(4)
*set var tmp8=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp9)
*cond(2) *cond(3) *cond(4)
*set var tmp9=-1
*set var nExcavation=operation(nExcavation+1)
*endif
*if(VolumeNumber==tmp10)
*cond(2) *cond(3) *cond(4)
*set var tmp10=-1
*set var nExcavation=operation(nExcavation+1)
*end if
*end elems
*end if
$$SUBMERGED_CALCULATION
*GenData(submerged_ID) *GenData(number_of_initialisation_steps_[-]__)
$$APPLY_OBJECTIVE_STRESS
*GenData(objective_stress_ID)
$$DEGREE_OF_FILLING
*GenData(degree_of_filling_[-]__)
$$NUMBER_OF_ACTIVE_ELEMENTS
*set var I=0
*loop elems onlyincond
*if((ElemsMat!=0))
*set var I=operation(I+1)
*end if
*end elems
*I
*##### DOUBLE-POINT FORMULATION #####
*if((GenData(doublepoint_ID,int))==1)
$$MAXIMUM_POROSITY
*GenData(Maximum_porosity_[-]__)
$$ERGUN_CONSTANTS
*GenData(Ergun_constants_-_a_[-]_and_b_[-]__)
$$PERMEABILITY_UPDATE
*if((GenData(permeabilityupdate_ID,int))==0)
*MessageBox **** INPUT ERROR: For the double-point formulation the permeability update has to be chosen. ****
*end if
*if((GenData(permeabilityupdate_ID,int))==1)
constant_permeability
$$INTRINSIC_PERMEABILITY
*GenData(intrinsic_permeability_[m2]__)
*end if
*if((GenData(permeabilityupdate_ID,int))==2)
Darcy_update
$$GRAIN_SIZE_DIAMETER
*GenData(Grain_size_diameter_[m]__)
*end if
*if((GenData(permeabilityupdate_ID,int))==3)
Ergun_update
$$GRAIN_SIZE_DIAMETER
*GenData(Grain_size_diameter_[m]__)
*end if
$$LIQUID_SURFACE
*GenData(Elevation_liquid_surface_for_K0_[m]__)
$$APPLY_STRAIN_SMOOTHING_LIQUID_TWOLAYERFORM
*GenData(Strain_smoothing_liquid__)
$$NO_TENSILE_STRESS_LIQUID_MP_WITH_LIQUID_STATUS
*GenData(No_tensile_stress_liquid__)
$$DETECT_FREE_SURFACE
*GenData(Detect_free_surface__)
*end if
*##### OUTPUT DATA #####
$$OUTPUT_NUMBER_OF_MATERIAL_POINTS
*if(strcasecmp(GenData(number_of_material_points__),"none")==0)
0
*else
*GenData(number_of_material_points__)
*end if
$$OUTPUT_MATERIAL_POINTS
*if(strcasecmp(GenData(number_of_material_points__),"none")==0)
-1
*end if
*if((GenData(number_of_material_points__,int))==1)
*GenData(material_point_#1__)
*end if
*if((GenData(number_of_material_points__,int))==2)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*end if
*if((GenData(number_of_material_points__,int))==3)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*end if
*if((GenData(number_of_material_points__,int))==4)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*end if
*if((GenData(number_of_material_points__,int))==5)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*end if
*if((GenData(number_of_material_points__,int))==6)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*GenData(material_point_#6__)
*end if
*if((GenData(number_of_material_points__,int))==7)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*GenData(material_point_#6__)
*GenData(material_point_#7__)
*end if
*if((GenData(number_of_material_points__,int))==8)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*GenData(material_point_#6__)
*GenData(material_point_#7__)
*GenData(material_point_#8__)
*end if
*if((GenData(number_of_material_points__,int))==9)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*GenData(material_point_#6__)
*GenData(material_point_#7__)
*GenData(material_point_#8__)
*GenData(material_point_#9__)
*end if
*if((GenData(number_of_material_points__,int))==10)
*GenData(material_point_#1__)
*GenData(material_point_#2__)
*GenData(material_point_#3__)
*GenData(material_point_#4__)
*GenData(material_point_#5__)
*GenData(material_point_#6__)
*GenData(material_point_#7__)
*GenData(material_point_#8__)
*GenData(material_point_#9__)
*GenData(material_point_#10__)
*end if
$$END
