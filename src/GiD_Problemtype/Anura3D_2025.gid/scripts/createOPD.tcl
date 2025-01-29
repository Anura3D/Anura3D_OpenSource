#      Write OPD file

proc Anura3D::WriteCalculationFile_OPD { stageNode projectPath projectName icount_stage} {             

    # RESULTS OUTPUT DATA  
    set filename [file join [file dir $projectPath] $projectName.A3D $projectName.OPD_stage${icount_stage}]
    if {[file exists $filename]} { file delete -force $filename }                       
    GiD_WriteCalculationFile init $filename  
    
    set xp {container[@n="results"]}
    set resultsNode [$stageNode selectNodes $xp]
    set xp {container[@n="output"]}
    set paraview_outputNode [$resultsNode selectNodes $xp]    
    set calc_steps_between_outputsNode \
        [$paraview_outputNode selectNodes {value[@n="calc_steps_between_outputs"]}]
    set calc_steps_between_outputs [$calc_steps_between_outputsNode @v]
    GiD_WriteCalculationFile puts {$$CALCULATION_STEPS_BETWEEN_OUTPUTS}
    GiD_WriteCalculationFile puts $calc_steps_between_outputs
    
    set output_varNode [$paraview_outputNode selectNodes {container[@n="output_variables"]}]
    
    set scalarNode [$output_varNode selectNodes {container[@n="scalar"]}]
    set generalNode [$scalarNode selectNodes {container[@n="general"]}]
     
    set mass_mixtureNode [$generalNode selectNodes {value[@n="mass_mixture"]}]
    set mass_mixture [$mass_mixtureNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$MASS_MIXTURE}
    GiD_WriteCalculationFile puts "$mass_mixture"
        
    set weight_mixtureNode [$generalNode selectNodes {value[@n="weight_mixture"]}]
    set weight_mixture [$weight_mixtureNode getAttribute "v"] 
    GiD_WriteCalculationFile puts {$$WEIGHT_MIXTURE}
    GiD_WriteCalculationFile puts "$weight_mixture"
    
    set porosityNode [$generalNode selectNodes {value[@n="porosity"]}]    
    set porosity [$porosityNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$POROSITY}
    GiD_WriteCalculationFile puts "$porosity"
    
    set mean_effective_stressNode [$generalNode selectNodes {value[@n="mean_effective_stress"]}]    
    set mean_effective_stress [$mean_effective_stressNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MEAN_EFFECTIVE_STRESS}
    GiD_WriteCalculationFile puts "$mean_effective_stress"
    
    set deviatoric_stressNode [$generalNode selectNodes {value[@n="deviatoric_stress"]}]    
    set deviatoric_stress [$deviatoric_stressNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$DEVIATORIC_STRESS}
    GiD_WriteCalculationFile puts "$deviatoric_stress"
    
    set deviatoric_strainNode [$generalNode selectNodes {value[@n="deviatoric_strain"]}]    
    set deviatoric_strain [$deviatoric_strainNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$DEVIATORIC_STRAIN}
    GiD_WriteCalculationFile puts "$deviatoric_strain"
    
    set degree_of_saturationNode [$generalNode selectNodes {value[@n="degree_of_saturation"]}]    
    set degree_of_saturation [$degree_of_saturationNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$DEGREE_OF_SATURATION}
    GiD_WriteCalculationFile puts "$degree_of_saturation"
    
    set state_variablesNode [$generalNode selectNodes {value[@n="state_variables"]}]   
    set state_variables [$state_variablesNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$STATE_VARIABLES}
    GiD_WriteCalculationFile puts "$state_variables"
    
    set dampingNode [$generalNode selectNodes {value[@n="damping"]}]   
    set damping [$dampingNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$DAMPING}
    GiD_WriteCalculationFile puts "$damping"
    
    set solidNode [$scalarNode selectNodes {container[@n="solid"]}]
    set mass_solidNode [$solidNode selectNodes {value[@n="mass_solid"]}]   
    set mass_solid [$mass_solidNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MASS_SOLID}
    GiD_WriteCalculationFile puts "$mass_solid"
    
    set weight_solidNode [$solidNode selectNodes {value[@n="weight_solid"]}]    
    set weight_solid [$weight_solidNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$WEIGHT_SOLID}
    GiD_WriteCalculationFile puts "$weight_solid"
    
    set volumetric_strain_solidNode [$solidNode selectNodes {value[@n="volumetric_strain_solid"]}]   
    set volumetric_strain_solid [$volumetric_strain_solidNode  getAttribute "v"]
    GiD_WriteCalculationFile puts {$$VOLUMETRIC_STRAIN_SOLID}
    GiD_WriteCalculationFile puts "$volumetric_strain_solid"            
    
    set liquidNode [$scalarNode selectNodes {container[@n="liquid"]}]
    set mass_liquidNode [$liquidNode selectNodes {value[@n="mass_liquid"]}]
    set mass_liquid [$mass_liquidNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MASS_LIQUID}
    GiD_WriteCalculationFile puts "$mass_liquid"
       
    set weight_liquidNode [$liquidNode selectNodes {value[@n="weight_liquid"]}]    
    set weight_liquid [$weight_liquidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$WEIGHT_LIQUID}
    GiD_WriteCalculationFile puts "$weight_liquid"
    
    set volumetric_strain_liquidNode [$liquidNode selectNodes {value[@n="volumetric_strain_liquid"]}]    
    set volumetric_strain_liquid [$volumetric_strain_liquidNode getAttribute "v"] 
    
    GiD_WriteCalculationFile puts {$$VOLUMETRIC_STRAIN_LIQUID}
    GiD_WriteCalculationFile puts "$volumetric_strain_liquid"
    
    set pressure_liquidNode [$liquidNode selectNodes {value[@n="pressure_liquid"]}]    
    set pressure_liquid [$pressure_liquidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$PRESSURE_LIQUID}
    GiD_WriteCalculationFile puts "$pressure_liquid"
    
    set gasNode [$scalarNode selectNodes {container[@n="gas"]}]
    
    set mass_gasNode [$gasNode selectNodes {value[@n="mass_gas"]}]    
    set mass_gas [$mass_gasNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MASS_GAS}
    GiD_WriteCalculationFile puts "$mass_gas"
    
    set weight_gasNode [$gasNode selectNodes {value[@n="weight_gas"]}]   
    set weight_gas [$weight_gasNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$WEIGHT_GAS}
    GiD_WriteCalculationFile puts "$weight_gas"
   
    set pressure_gasNode [$gasNode selectNodes {value[@n="pressure_gas"]}]    
    set pressure_gas [$pressure_gasNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$PRESSURE_GAS}
    GiD_WriteCalculationFile puts "$pressure_gas" 
    
    set volumetric_strain_gasNode [$gasNode selectNodes {value[@n="volumetric_strain_gas"]}]   
    set volumetric_strain_gas [$volumetric_strain_gasNode getAttribute "v"] 
    GiD_WriteCalculationFile puts {$$VOLUMETRIC_STRAIN_GAS}
    GiD_WriteCalculationFile puts "$volumetric_strain_gas"      

    # VECTOR    
    set vectorNode [$output_varNode selectNodes {container[@n="vector"]}]
    set generalNode [$vectorNode selectNodes {container[@n="general"]}]
    set global_posNode [$generalNode selectNodes {value[@n="global_pos"]}]
    set global_pos [$global_posNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$GLOBAL_POSITION}    
    GiD_WriteCalculationFile puts "$global_pos"
    
    set local_posNode [$generalNode selectNodes {value[@n="local_pos"]}]    
    set local_pos [$local_posNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$LOCAL_POSITION}
    GiD_WriteCalculationFile puts "$local_pos"
    
    set body_force_mixtureNode [$generalNode selectNodes {value[@n="body_force_mixture"]}]   
    set body_force_mixture [$body_force_mixtureNode getAttribute "v"] 
    GiD_WriteCalculationFile puts {$$BODY_FORCE_MIXTURE}
    GiD_WriteCalculationFile puts "$body_force_mixture"

    set solidNode [$vectorNode selectNodes {container[@n="solid"]}]
    set body_force_solidNode [$solidNode selectNodes {value[@n="body_force_solid"]}]   
    set body_force_solid [$body_force_solidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$BODY_FORCE_SOLID}
    GiD_WriteCalculationFile puts "$body_force_solid"
    
    set external_force_solidNode [$solidNode selectNodes {value[@n="external_force_solid"]}]    
    set external_force_solid [$external_force_solidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$EXTERNAL_FORCE_SOLID}
    GiD_WriteCalculationFile puts "$external_force_solid"
        
    set acceleration_solidNode [$solidNode selectNodes {value[@n="acceleration_solid"]}]    
    set acceleration_solid [$acceleration_solidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$ACCELERATION_SOLID}
    GiD_WriteCalculationFile puts "$acceleration_solid"
      
    set velocity_solidNode [$solidNode selectNodes {value[@n="velocity_solid"]}]    
    set velocity_solid [$velocity_solidNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$VELOCITY_SOLID}
    GiD_WriteCalculationFile puts "$velocity_solid"
    
    set disp_solidNode [$solidNode selectNodes {value[@n="disp_solid"]}]    
    set disp_solid [$disp_solidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$DISPLACEMENT_SOLID}
    GiD_WriteCalculationFile puts "$disp_solid"

    set liquidNode [$vectorNode selectNodes {container[@n="liquid"]}]  
    set body_force_liquidNode [$liquidNode selectNodes {value[@n="body_force_liquid"]}]    
    set body_force_liquid [$body_force_liquidNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$BODY_FORCE_LIQUID}
    GiD_WriteCalculationFile puts "$body_force_liquid"
        
    set external_force_liquidNode [$liquidNode selectNodes {value[@n="external_force_liquid"]}]   
    set external_force_liquid [$external_force_liquidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$EXTERNAL_FORCE_LIQUID}
    GiD_WriteCalculationFile puts "$external_force_liquid" 
        
    set velocity_liquidNode [$liquidNode selectNodes {value[@n="velocity_liquid"]}]    
    set velocity_liquid [$velocity_liquidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$VELOCITY_LIQUID}
    GiD_WriteCalculationFile puts "$velocity_liquid"
    
    set disp_liquidNode [$liquidNode selectNodes {value[@n="disp_liquid"]}]    
    set disp_liquid [$disp_liquidNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$DISPLACEMENT_LIQUID}
    GiD_WriteCalculationFile puts "$disp_liquid"
        
    set gasNode [$vectorNode selectNodes {container[@n="gas"]}]  
    set body_force_gasNode [$gasNode selectNodes {value[@n="body_force_gas"]}]    
    set body_force_gas [$body_force_gasNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$BODY_FORCE_GAS}
    GiD_WriteCalculationFile puts "$body_force_gas"
        
    set external_force_gasNode [$gasNode selectNodes {value[@n="external_force_gas"]}]   
    set external_force_gas [$external_force_gasNode getAttribute "v"]
    
    GiD_WriteCalculationFile puts {$$EXTERNAL_FORCE_GAS}
    GiD_WriteCalculationFile puts "$external_force_gas" 
        
    set velocity_gasNode [$gasNode selectNodes {value[@n="velocity_gas"]}]    
    set velocity_gas [$velocity_gasNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$VELOCITY_GAS}
    GiD_WriteCalculationFile puts "$velocity_gas"
    
    set disp_gasNode [$gasNode selectNodes {value[@n="disp_gas"]}]    
    set disp_gas [$disp_gasNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$DISPLACEMENT_GAS}
    GiD_WriteCalculationFile puts "$disp_gas"
    
        # TENSOR    
    set tensorNode [$output_varNode selectNodes {container[@n="tensor"]}]
    set total_stressNode [$tensorNode selectNodes {value[@n="total_stress"]}]          
    set total_stress [$total_stressNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$TOTAL_STRESS}
    GiD_WriteCalculationFile puts "$total_stress"
    
    set effective_stressNode [$tensorNode selectNodes {value[@n="effective_stress"]}]        
    set effective_stress [$effective_stressNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$EFFECTIVE_STRESS}
    GiD_WriteCalculationFile puts "$effective_stress"
    
    set strainNode [$tensorNode selectNodes {value[@n="strain"]}]         
    set strain [$strainNode getAttribute "v"]
    GiD_WriteCalculationFile puts {$$STRAIN}
    GiD_WriteCalculationFile puts "$strain"
        
    # Elements for material point output
    set num_elem 0  
    set ov_type "element"    
    set xp_res [format_xpath {container[@n="results"]/container[@n="elements_for_mat_point_output"]/condition[@n="select_elements_for_mat_point_output"]/group[@ov=%s]} $ov_type]    
    foreach gNode [$stageNode selectNodes $xp_res] {
            set gName [get_domnode_attribute $gNode n]
            set gEntities_num [GiD_EntitiesGroups get $gName $ov_type -count]
            set num_elem [expr $gEntities_num + $num_elem] 
    }
		
	set xp_path {container[@n="results"]/container[@n="elements_for_mat_point_output"]}
    set xp [$stageNode selectNodes $xp_path]
	set num_output_MP [$xp selectNodes {string(value[@n="number_output_mat_point"]/@v)}]
	
    set format ""
    foreach gNode [$stageNode selectNodes $xp_res] {
            set n [$gNode @n]
            dict set format $n "%d \n"
    }
        
        if { [dict size $format] } {
            set err [catch {
                GiD_WriteCalculationFile elements -return $format} numelems]       
            if { $err } {          
                snit_messageBox -parent .gid \
                    -message [= "Error when writing elements for material point output (%s)" $numelems]             
                return "isnotok"                        
            } else {    
                if { $icount_stage == "0"} {
				GiD_WriteCalculationFile puts {$$NUMBER_OF_OUTPUT_MP_PER_ELEMENT}
                GiD_WriteCalculationFile puts $num_output_MP	
                GiD_WriteCalculationFile puts {$$NUMBER_OF_ELEMENTS_FOR_MP_OUTPUT}
                GiD_WriteCalculationFile puts $num_elem
                GiD_WriteCalculationFile puts {$$ELEMENTS_FOR_MP_OUTPUT}        
                GiD_WriteCalculationFile elements $format
            }
        }                				            
    }        

	set time_between_outputs [$xp selectNodes {string(value[@n="time_between_outputs"]/@v)}]   
    GiD_WriteCalculationFile puts {$$TIME_BETWEEN_MP_OUTPUTS}
    GiD_WriteCalculationFile puts $time_between_outputs
			
    GiD_WriteCalculationFile puts {$$END}
    GiD_WriteCalculationFile end   
}
