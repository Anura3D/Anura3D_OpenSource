proc write_material_properties { dimension_type root } {
    set xp [format_xpath {container[@n="materials"]/blockdata}]
    set list [$root selectNodes $xp]
    set list_len [llength $list]
    GiD_WriteCalculationFile puts {$$NUMBER_OF_MATERIALS}
    GiD_WriteCalculationFile puts $list_len
    set int 1
    foreach gNode $list {

        set type [$gNode getAttribute "name"]
        set MATERIAL_ID [find_material_id $type $root]
        GiD_WriteCalculationFile puts {$$MATERIAL_INDEX}
        GiD_WriteCalculationFile puts $MATERIAL_ID
        GiD_WriteCalculationFile puts {$$MATERIAL_NAME}
        GiD_WriteCalculationFile puts $type
        set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="material_type_"]}]]
        set typename [$node getAttribute "v"]

        if {$typename == "Dry material"} {
            set type "dry_material"
        } elseif {$typename == "Saturated material-drained"} {
            set type "saturated_material_drained"
        } elseif {$typename == "Saturated material-undrained effective stress"} {
            set type "saturated_material_undrained_effective"
        } elseif {$typename == "Saturated material-undrained total stress"} {
            set type "saturated_material_undrained_total"
        } elseif {$typename == "Saturated material-fully coupled"} {
            set type "saturated_material_coupled"
        } elseif {$typename == "Unsaturated material-2-phase with suction effect"} {
            set type "unsaturated_material_2phase_suction"
        } elseif {$typename == "Unsaturated material-3-phase fully coupled"} {
            set type "unsaturated_material_2phase_suction"
        } elseif {$typename == "Liquid"} {
            set type "liquid"
        }

        GiD_WriteCalculationFile puts {$$MATERIAL_TYPE}
        GiD_WriteCalculationFile puts $type

        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="initial_porosity_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$POROSITY_SOLID}
            GiD_WriteCalculationFile puts $type
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_solid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_SOLID }
            GiD_WriteCalculationFile puts $type
        }
        if {$typename != "Dry material"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_LIQUID }
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="density_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DENSITY_GAS}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Saturated material-fully coupled"||$typename == "Unsaturated material-2-phase with suction effect"||$typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="intrinsic_permeability_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY_LIQUID}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="intrinsic_permeability_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$INTRINSIC_PERMEABILITY_GAS}
            GiD_WriteCalculationFile puts $type
        }
        if {($typename == "Liquid") || ($typename == "Saturated material-fully coupled") || ($typename == "Unsaturated material-2-phase with suction effect") || ($typename == "Unsaturated material-3-phase fully coupled")} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="bulk__modulus_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$BULK_MODULUS_LIQUID}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="bulk__modulus_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$BULK_MODULUS_GAS}
            GiD_WriteCalculationFile puts $type
        }
        if {($typename == "Liquid") || ($typename == "Saturated material-fully coupled") || ($typename == "Unsaturated material-2-phase with suction effect") || ($typename == "Unsaturated material-3-phase fully coupled") } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="dinamic_viscosity_liquid_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DYNAMIC_VISCOSITY_LIQUID}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Unsaturated material-3-phase fully coupled" } {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="dinamic_viscosity_gas_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$DYNAMIC_VISCOSITY_GAS}
            GiD_WriteCalculationFile puts $type
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="elastic_swelling_index_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$SWELLING_INDEX}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="K0-value_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$K0_VALUE_SOLID}
            GiD_WriteCalculationFile puts $type
        }
        if {$typename == "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="liquid_cavitation_"]}]]
            set type [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$LIQUID_CAVITATION}
            GiD_WriteCalculationFile puts $type
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="detect_free_surface_liquid_"]}]]
            set type [$node getAttribute "v"]
            set node [$gNode selectNodes [format_xpath {container[@n="_basic"]/value[@n="free_surface_factor_"]}]]
            set factor [$node getAttribute "v"]

            GiD_WriteCalculationFile puts {$$APPLY_DETECT_LIQUID_SURFACE}
            if {$type == "No"} {
                GiD_WriteCalculationFile puts [= "0 %s" $factor]
            } elseif {$type == "Yes"} {
                GiD_WriteCalculationFile puts [= "1 %s" $factor]
            }
        }

        if {$typename != "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_material_model_solid_"]}]]
            set typemodel [$node getAttribute "v"]
            if {$typemodel == "Rigid body"} {
                set model "rigid_body"
            } elseif {$typemodel == "Linear Elasticity"} {
                set model "linear_elasticity"
            } elseif {$typemodel == "Mohr-Coulomb"} {
                set model "mohr_coulomb"
            } elseif {$typemodel == "External Material Model"} {
                set model "external_soil_model"
            }

            GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_SOLID}
            GiD_WriteCalculationFile puts $model

            if {$typemodel == "Rigid body"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="x-constr"]}]]
                set type [$node getAttribute "v"]
                if {$type == "Yes"} {
                    set flag 1
                } else {
                    set flag 0
                }
                GiD_WriteCalculationFile puts {$$constraint_XDISPLACEMENT}
                GiD_WriteCalculationFile puts $flag
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="y-constr"]}]]
                set type [$node getAttribute "v"]
                if {$type == "Yes"} {
                    set flag 1
                } else {
                    set flag 0
                }
                GiD_WriteCalculationFile puts {$$constraint_YDISPLACEMENT}
                GiD_WriteCalculationFile puts $flag
                if {$dimension_type == "3D" || $dimension_type == "3D:Axissymmetric"} {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="z-constr"]}]]
                    set type [$node getAttribute "v"]
                    if {$type == "Yes"} {
                        set flag 1
                    } else {
                        set flag 0
                    }
                    GiD_WriteCalculationFile puts {$$constraint_ZDISPLACEMENT}
                    GiD_WriteCalculationFile puts $flag
                }
            } elseif {$typemodel == "Linear Elasticity"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_Young_modulus_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$YOUNG_MODULUS}
                GiD_WriteCalculationFile puts $type
                if {$typename != "Saturated material-undrained total stress"} {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
                    set type [$node getAttribute "v"]
                    GiD_WriteCalculationFile puts {$$POISSON_RATIO}
                    GiD_WriteCalculationFile puts $type
                } else {
                    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
                    set type [$node getAttribute "v"]
                    GiD_WriteCalculationFile puts {$$POISSON_RATIO}
                    GiD_WriteCalculationFile puts $type
                }
            } elseif {$typemodel == "Mohr-Coulomb"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_Young_modulus_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$YOUNG_MODULUS}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_poisson_ratio_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$POISSON_RATIO}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_friction_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$FRICTION_ANGLE}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_cohesion_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$COHESION}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_dilatancy_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$DILATANCY_ANGLE}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="eff_tensile_strength_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$TENSILE_STRENGTH}
                GiD_WriteCalculationFile puts $type

            } elseif {$typemodel == "External Material Model"} {
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_model_dll_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_DLL}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_model_dll_dim_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$UMAT_DIMENSION}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_01_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_01}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_02_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_02}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_03_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_03}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_04_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_04}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_05_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_05}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_06_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_06}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_07_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_07}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_08_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_08}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_09_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_09}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_10_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_10}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_11_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_11}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_12_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_12}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_13_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_13}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_14_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_14}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_15_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_15}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_16_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_16}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_17_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_17}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_18_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_18}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_19_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_19}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_20_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_20}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_21_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_21}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_22_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_22}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_23_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_23}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_24_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_24}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_25_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_25}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_26_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_26}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_27_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_27}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_28_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_28}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_29_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_29}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_30_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_30}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_31_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_31}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_32_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_32}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_33_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_33}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_34_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_34}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_35_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_35}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_36_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_36}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_37_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_37}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_38_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_38}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_39_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_39}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_40_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_40}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_41_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_41}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_42_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_42}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_43_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_43}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_44_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_44}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_45_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_45}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_46_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_46}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_47_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_47}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_48_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_48}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_49_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_49}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="material_parameter_solid_50_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_50}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_01_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_01}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_02_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_02}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_03_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_03}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_04_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_04}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_05_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_05}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_06_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_06}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_07_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_07}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_08_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_08}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_09_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_09}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_10_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_10}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_11_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_11}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_12_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_12}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_13_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_13}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_14_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_14}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_15_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_15}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_16_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_16}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_17_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_17}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_18_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_18}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_19_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_19}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_20_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_20}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_21_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_21}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_22_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_22}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_23_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_23}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_24_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_24}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_25_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_25}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_26_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_26}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_27_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_27}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_28_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_28}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_29_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_29}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_30_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_30}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_31_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_31}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_32_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_32}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_33_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_33}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_34_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_34}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_35_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_35}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_36_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_36}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_37_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_37}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_38_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_38}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_39_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_39}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_40_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_40}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_41_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_41}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_42_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_42}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_43_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_43}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_44_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_44}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_45_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_45}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_46_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_46}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_47_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_47}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_48_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_48}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_49_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_49}
                GiD_WriteCalculationFile puts $type

                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="initial_state_variable_solid_50_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$INITIAL_STATE_VARIABLE_SOLID_50}
                GiD_WriteCalculationFile puts $type

            }
        }
        if {$typename == "Liquid"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_material_model_liquid_"]}]]
            set typemodel [$node getAttribute "v"]
            GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_LIQUID}
            if {$typemodel == "Newtonian"} {
                GiD_WriteCalculationFile puts  "newtonian_liquid"
            } elseif {$typemodel == "Bingham Fluid"} {
                GiD_WriteCalculationFile puts  "bingham_liquid"
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="bingham_yield_stress_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_YIELD_STRESS}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_young_modulus_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_YOUNG_MODULUS}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_poisson_ratio_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$BINGHAM_POISSON_RATIO}
                GiD_WriteCalculationFile puts $type
            } elseif {$typemodel == "Frictional Fluid"} {
                GiD_WriteCalculationFile puts  "frictional_liquid"
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="fluid_friction_angle_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_FRICTION_ANGLE}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_young_modulus_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_YOUNG_MODULUS}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="elastic_poisson_ratio_l_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$LIQUID_POISSON_RATIO}
                GiD_WriteCalculationFile puts $type
            }
        }

        if {$typename == "Unsaturated material-2-phase with suction effect" || $typename == "Unsaturated material-3-phase fully coupled"} {
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_retention_curve_"]}]]
            set type [$node getAttribute "v"]
            if {$type == "Linear"} {
                GiD_WriteCalculationFile puts {$$WATER_RETENTION_CURVE}
                GiD_WriteCalculationFile puts {linear}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_rc_linear"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$av}
                GiD_WriteCalculationFile puts $type
            } elseif {$type == "Van Genuchten"} {
                GiD_WriteCalculationFile puts {$$WATER_RETENTION_CURVE}
                GiD_WriteCalculationFile puts {van_genuchten}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_min_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smin}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_max_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smax}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_ref_press_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$P0}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="slrc_lambda_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$L}
                GiD_WriteCalculationFile puts $type
            }
            set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="_unsat_hydraulic_cond_"]}]]
            set type [$node getAttribute "v"]

            if {$type == "Constant"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {constant}
            } elseif {$type == "Hillel"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {hillel}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_r_exp"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$r}
                GiD_WriteCalculationFile puts $type
            } elseif {$type == "Mualem"} {
                GiD_WriteCalculationFile puts {$$HYDR_CONDUCTIVITY_CURVE}
                GiD_WriteCalculationFile puts {mualem}
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_min_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smin}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_max_deg_sat_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$Smax}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_ref_press_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$P0}
                GiD_WriteCalculationFile puts $type
                set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="hycon_lambda_"]}]]
                set type [$node getAttribute "v"]
                GiD_WriteCalculationFile puts {$$L}
                GiD_WriteCalculationFile puts $type
            }
        }
    }
}