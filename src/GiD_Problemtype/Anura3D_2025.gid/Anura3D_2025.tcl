set err [catch { package present Tk }]
if { !$err } {
    set ::has_tk 1
} else {
    set ::has_tk 0
}

if { [info commands ::GidUtils::IsTkDisabled] ne "" && [GidUtils::IsTkDisabled] } {
    set ::has_tk 0
}

package require gid_cross_platform

proc InitGIDProject { dir } {

    global _dir
    global ProgramName VersionNumber Priv
    
    set _dir $dir
    set anura_version "Anura3D v2025"
    set GiDVersionRequired "14.0"

    Anura3D::SetDir $dir ;#store to use it later
    Anura3D::LoadScripts
    GidUtils::OpenWindow CUSTOMLIB
    
    dom parse [tDOM::xmlReadFile [file join $dir Anura3D_2025.xml]] doc
    set ProgramName [$doc selectNodes string(Infoproblemtype/Program/Name)]
    set VersionNumber [$doc selectNodes string(Infoproblemtype/Program/Version)]
       
    set Priv(problemtypedir) $dir

    # disclaimer at start-up
    set self_close 0
    GidUtils::Splash [file join $_dir images about_anura3d.png] .splash $self_close [list $anura_version 445 10]
    
    # check minimal GiD version
    set message_required "This $anura_version interface is developed for GiD $GiDVersionRequired or later. \n \n It is advised to update your GiD software."
    set title_required "Anura3D - Required GiD Version"
    if { [GidUtils::VersionCmp $GiDVersionRequired] < 0 } { tk_messageBox -title $title_required -message $message_required -icon warning -type ok }
   
    # Add ANURA3D menu
    GiDMenu::Create "Anura3D" "PRE" 5 =
    GiDMenu::InsertOption "Anura3D" [list "Calculate"] 0 PRE "Anura3D::Calculate local" "" "[file join $dir images calculate.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "Stages manager"] 1 PRE "Anura3D::StagesManager" "" "[file join $dir images configure.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "View process info"] 2 PRE "Anura3D::ViewProcessInfo" "" "[file join $dir images info.png]" insert =  
    GiDMenu::InsertOption "Anura3D" [list "---"] 3 PRE "" "" "" insert =    
    GiDMenu::InsertOption "Anura3D" [list "Tutorial Manual..."] 4 PRE "Anura3D::Tutorial" "" "[file join $dir images tutorial.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "Scientific Manual..."] 5 PRE "Anura3D::Scientific" "" "[file join $dir images scientific.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "Verification Manual..."] 6 PRE "Anura3D::Verification" "" "[file join $dir images verification.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "Disclaimer..."] 7 PRE "Anura3D::Disclaimer" "" "[file join $dir images icon_anura3d.png]" insert =
    GiDMenu::InsertOption "Anura3D" [list "About..."] 8 PRE "Anura3D::About" "" "[file join $dir images icon_anura3d.png]" insert =

    # simplify and adaptd HELP menu   
    GiDMenu::RemoveOption "Help" [list "Tutorials"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "What is new"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "FAQ"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "Register problem type"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "Register from file"] "PRE" _

    # remove options in DATA menu
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Local axes" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Problem data" ""
    
    # remove CALCULATE menu
    GiDMenu::Delete "Calculate" PRE

    # apply menu update
    Anura3D::Bitmaps $dir

    GiDMenu::UpdateMenus

    # open window tree
    GidUtils::OpenWindow CUSTOMLIB

}

proc EndGIDProject { } {

}

namespace eval Anura3D {
}

namespace eval ANURA3D_2025 {
    variable problemtype_dir
}

proc Anura3D::SetDir { dir } {
    variable problemtype_dir
    set problemtype_dir $dir
}

proc Anura3D::GetDir { } {
    variable problemtype_dir
    return $problemtype_dir
}

proc Anura3D::Tutorial { } {

    global _dir
    set TestDoc [file join $_dir doc "TutorialManual_2025.pdf"]
    eval exec [auto_execok start] \"\" [list $TestDoc]

}

proc Anura3D::Scientific { } {

    global _dir
    set TestDoc [file join $_dir doc "ScientificManual_2022.pdf"]
    eval exec [auto_execok start] \"\" [list $TestDoc]

}

proc Anura3D::Verification { } {

    global _dir
    set TestDoc [file join $_dir doc "VerificationManual_2021.pdf"]
    eval exec [auto_execok start] \"\" [list $TestDoc]

}

proc Anura3D::ViewProcessInfo {} {
    global GidProcWin RunProcInfo
   
    if { ![info exists RunProcInfo] } {
        set RunProcInfo "run"
    }

    if { ![info exists ::GidProcWin(w)] || ![winfo exists $::GidProcWin(w).listbox#1] } {
        set wbase .gid
        set w ""
    } else {
        set wbase $::GidProcWin(w)
        set w $::GidProcWin(w).listbox#1
    }
   
    set project_name [GiD_Info Project ModelName]
    if { [file extension $project_name] == ".gid" } {
        set project_name [file root $project_name]
    }
    set basename [file tail $project_name]    
    set filename [file join $project_name.A3D $basename.PROC]           
    
    PWViewOutputWin $filename current "" ""        
} 

proc Anura3D::About { } {

    snit_messageBox -parent .gid -title "Anura3D - About" \
        -message "Anura3D Version 2025\n\nFor more information about Anura3D, check the website\n\nhttp://www.anura3D.com" -icon info

}

proc Anura3D::Disclaimer { } {

    set answer [snit_messageBox -parent .gid -title "Anura3D - Disclaimer" -message "Copyright (C) 2020  Members of the Anura3D MPM Research Community

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>

Anura3D MPM Research Community
E-mail:        info@anura3D.com
Web:         www.anura3D.com" -icon info -type ok]

}

proc Anura3D::SelectBatFile { } {     
    
    set ret ""
    switch $::tcl_platform(platform) {
        "windows" {     
            set ret [list Anura3D.bat]          
        }
        "default" { 
            WarnWin [= {This version cannot run the Anura 3D calculation}]
            return ""
        }
    }    
    return $ret
}

proc Anura3D::StageToCalculate { stagename stagestocalcList } {
    set stagetocalculate 0   
    foreach i [lindex $stagestocalcList 0] {
        if { $i eq $stagename } { 
            set stagetocalculate 1 
            break
        }
    }   
    return $stagetocalculate
}

proc Anura3D::Calculate { args } {
    
    global Priv
    
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]  
    
    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    } 
    
    set dir $Priv(problemtypedir)
    set projectPath [GiD_Info Project ModelName]
    set projectfolder "$projectPath.gid"
    set projectName [file tail $projectPath]     

    if { $::has_tk } {
        if { $projectName eq "UNNAMED" } {
            snit_messageBox -parent .gid \
                -message [= "Before calculating, a project title is needed. Save project to get it"]
            return   
        }            
        if { ![llength [GiD_Geometry list point 1:]] } {
            snit_messageBox -parent .gid -message \
                [= "No mesh generated. No geometry present"]
            return       
        }     
    } else {
        if { $projectName eq "UNNAMED" } {
            error[= "Before calculating, a project title is needed. Save project to get it"]
            return   
        }            
        if { ![llength [GiD_Geometry list point 1:]] } {
            error [= "No mesh generated. No geometry present"]
            return       
        }
    }       
    
    set title [= "Calculate"]
    set message [= "Before creating the Anura3D input files, make sure that: \n \n \
            - Materials are assigned \n  - Boundary conditions are defined \n  - Material points are specified \n  - Calculation parameters are defined \n \
            - Mesh is (re-)generated.\n\nSave the project files and make sure that the project name\ndoes NOT contain any spaces!"]
    
    set answer [snit_messageBox -title $title -parent .gid -message $message \
        -icon warning -type okcancel]

    switch -- $answer {
        cancel {
            return  
        }
        ok {      
            set A3D_folder [file join [file dir $projectPath] $projectName.A3D]
            if { $args eq "local" } {
                if { [file exists $A3D_folder] } {
                    file delete -force $A3D_folder
                }
            } else {
                set calculateallstages 1
                set xp {container[@n='stages']/blockdata[@n='stage']} 
                foreach stageNode [$root selectNodes $xp]  {
                    set stagetocalculate [Anura3D::StageToCalculate [$stageNode @name] $args]
                    if { !$stagetocalculate } {
                        set calculateallstages 0
                    } 
                }
                if { $calculateallstages } {
                    if { [file exists $A3D_folder] } {
                        file delete -force $A3D_folder
                    }
                }
            }               
            if { ![file exists $A3D_folder] } { file mkdir $A3D_folder }           
            
            file copy -force [file join $dir dll zlib.dll] [file join $A3D_folder zlib.dll]
            file copy -force [file join $dir dll zlib.dll] [file join $A3D_folder libiomp5md.dll]
            file copy -force [file join $dir exec Anura3D_2025.exe] [file join $A3D_folder Anura3D_2025.exe]
            
            # For each stage              
            set xp {container[@n='stages']/blockdata[@n='stage']} 
            set j 0
            foreach stageNode [$root selectNodes $xp]  {                            
                set last_stage $j
                incr j                                
            }       
            
            set k 0
            foreach stageNode [$root selectNodes $xp]  {
                if {$k == $last_stage} {
                    set xp_tot_num_steps {container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/}
                    append xp_tot_num_steps {value[@n="number_of_calculation_steps"]}                            
                    set tot_num_stepsNode [$stageNode selectNodes $xp_tot_num_steps]                           
                    set tot_num_steps [$tot_num_stepsNode @v]                             
                    break        
                }
                incr k
            }        
                        
            set icount_stage 0
            foreach stageNode [$root selectNodes $xp]  {                 
                set stagetocalculate 0                
                if { $args eq "local" } {
                    # General calculation of all the stages
                    set stagetocalculate 1
                } else {
                    # Stages manager calculation
                    set stagetocalculate [Anura3D::StageToCalculate [$stageNode @name] $args]
                }                    
                if { $stagetocalculate } {  

                    Anura3D::ViewProcessInfo                    
                    
                    # Write PSF_$num files 
                    Anura3D::WritePSFFiles $stageNode $projectPath $projectName 
                    
                    # Write HHBF_$num files 
                    Anura3D::WriteHHBFFile $stageNode $projectPath $projectName
                    
                    # Write PVF file
                    Anura3D::WritePVFFile $stageNode $projectPath $projectName
                    
                    # Write CPS file
                    if { $icount_stage == "0"} {                       
                        set cps_extension "CPS_001"
                        set prev_stage_steps "0"
                        set total_time "0.0" 
                        set out_num_material_points "0"
                        set out_material_points "0"
                        set i 0
                        set layer_path {string(container[@n="General_data"]/value[@n="NLAYERS"]/@v)}
                        set layer_type [$root selectNodes $layer_path]
                        set xp_dim {container[@n="General_data"]/value[@n="NDIM"]}
                        set dim_typeNode [$root selectNodes $xp_dim]
                        set dim_type [$dim_typeNode @v]
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
                        foreach gNode [$stageNode selectNodes $xp] {                       
                            set l_material [$gNode selectNodes {string(value[@n="material"]/@v)}]
                            set list_group [$gNode @n]
                            set elements_id [GiD_EntitiesGroups get $list_group elements]
                            set num_elems [objarray length $elements_id]
                            set i [expr $i + $num_elems]
                        }
                        set number_of_active_elements $i
                    } else {
                        set ii 0
                        set stage_xp {container[@n='stages']/blockdata[@n='stage']} 
                        foreach istageNode [$root selectNodes $stage_xp] {
                            if {$ii == [expr {$icount_stage - 1}]} {
							    set cps_extension [ReadGeneralInfo cps_extension_after_calc $istageNode]
                                set prev_stage_steps [ReadGeneralInfo prev_stage_steps $istageNode]
                                set total_time [ReadGeneralInfo total_time $istageNode]
                                set number_of_active_elements [ReadGeneralInfo number_of_active_elements $istageNode]
								set out_num_material_points [ReadGeneralInfo out_num_material_points $istageNode]
								set out_material_points [ReadGeneralInfo out_material_points $istageNode]                                
                                break        
                            }
                        incr ii
                        }
                    }   
                                        
                    # REMOVE FILES FROM THE CURRENT AND SUBSEQUENT STAGES IF THEY EXIST
                    set first_step_stage [expr {$prev_stage_steps + 1}]
                    if { $first_step_stage < 10 } {
                        set file_ext "_00$first_step_stage"
						set file_ext2 "00$first_step_stage"
                    } else {
                        set file_ext "_0$first_step_stage"
						set file_ext2 "0$first_step_stage"
                    }                    

					for {set iii $first_step_stage} {$iii <= $tot_num_steps} {incr iii} {
                        if { $iii < 10 } {
							set ext "_00$iii"
                            } else {
                            set ext "_0$iii"
						}
						set BRF_file [file join $A3D_folder $projectName.BRF$ext]
                        if { [file exists $BRF_file] } {
							file delete $BRF_file
						}
						set GIP_file [file join $A3D_folder $projectName.GIP$ext]
                        if { [file exists $GIP_file] } {
						file delete $GIP_file	
						}
						set ENG_file [file join $A3D_folder $projectName$ext.ENG]
                        if { [file exists $ENG_file] } {
						file delete $ENG_file
						}
						set INF_file [file join $A3D_folder $projectName$ext.INF]
                        if { [file exists $INF_file] } {
						file delete $INF_file
						}
						set MLG_file [file join $A3D_folder $projectName$ext.MLG]
                        if { [file exists $MLG_file] } {
						file delete $MLG_file
						}
						set CPS_file [file join $A3D_folder $projectName.CPS$ext]
                        if { [file exists $CPS_file] } {
						file delete $CPS_file
						}
						set meshdata "_MeshData"
						set VTK_meshdata_file [file join $A3D_folder $projectName$meshdata$ext.vtk]
                        if { [file exists $VTK_meshdata_file] } {
						file delete $VTK_meshdata_file
						}
						set scalar "_MPScalar"
						set VTK_scalar_file [file join $A3D_folder $projectName$scalar$ext.vtk]
						if { [file exists $VTK_scalar_file] } {
						file delete $VTK_scalar_file
						}
						set tensor "_MPTensor"
						set VTK_tensor_file [file join $A3D_folder $projectName$tensor$ext.vtk]
                        if { [file exists $VTK_tensor_file] } {
						file delete $VTK_tensor_file
						}
						set vector "_MPVector"
						set VTK_vector_file [file join $A3D_folder $projectName$vector$ext.vtk]
                        if { [file exists $VTK_vector_file] } {
						file delete $VTK_vector_file
						}
					}	                  
                                        
                    # REMOVE GiD POST-PROCESS FILES FROM THE CURRENT AND SUBSEQUENT STAGES IF THEY EXIST
                    set POST_bin_file [file join $A3D_folder $projectName$file_ext2.POST.bin]
                    if { [file exists $POST_bin_file] } {
                        for {set iii $first_step_stage} {$iii <= $tot_num_steps} {incr iii} {
                            if { $iii < 10 } {
                                set ext "00$iii"
                            } else {
                                set ext "0$iii"
                            }
                        set POST_bin_file [file join $A3D_folder $projectName$ext.POST.bin]
                        file delete $POST_bin_file
                        }
                    }
                                        
                    # Write CPS file
                    set CPS_name [file join $projectfolder "$projectName-1.dat"]    
                    if { [file exists $CPS_name] } {
                        file delete $CPS_name
                    }              
                    Anura3D::WriteCalculationFile_CPS $CPS_name $stageNode $icount_stage $total_time $number_of_active_elements $out_num_material_points $out_material_points
                    
                    # Write GOM file               
                    set GOM_name [file join $projectfolder "$projectName-2.dat"]                          
                    if { [file exists $GOM_name] } {
                        file delete $GOM_name
                    }
                    Anura3D::WriteCalculationFile_GOM $GOM_name $stageNode $projectPath $projectName
                    set gom_extension "GOM_stage$icount_stage"
                    
                    # Write OPD_stage$number file with results output data
					Anura3D::WriteCalculationFile_OPD $stageNode $projectPath $projectName $icount_stage   
					
                    if {![file exists [file join $Priv(problemtypedir) exec Anura3D_2025.exe]]} {
                        snit_messageBox -parent .gid \
                            -message "This distribution does not contain the executable (.exe)"
                        return
                    }   
                    
                    set batfilename [Anura3D::SelectBatFile] 
                    set batchfile [file join $projectfolder $batfilename]
                    set err [catch {Anura3D::anura3D_WriteBatchFile $batchfile $icount_stage \
                                $cps_extension $gom_extension} errstring]                
                    if { $err } {
                        snit_messageBox -parent .gid -message "$errstring"
                        return
                    }                                   
                    exec $batchfile $projectName [file nativename $projectfolder] \
                        [file nativename $dir] [file nativename $A3D_folder] 
                        
                    # Call .exe executable                    
                    set pid [exec [file join $A3D_folder Anura3D_2025.exe] [file join $A3D_folder $projectName] PAUSE &]                      
                    
                    while 1 {                      
                        if { ![isalive $pid] } { break }
                        after 200
                        update
                    }
                    
                    set POST_file [file join $A3D_folder $projectName.POST.lst]                          
                    if { [file exists $POST_file] } {
                        file copy -force [file join $A3D_folder $projectName.POST.lst] [file join $projectfolder $projectName.POST.lst]
                    }
	
                    set err_file [file join $A3D_folder $projectName.err]
                    
                    if { [file exists $err_file] } {
                        set fileSize [file size $err_file]
                        if { $fileSize > 0 } {
                            UpdateGeneralInfo calculation_state "wrongcalc" $stageNode
                            snit_messageBox -parent .gid -message "$errstring"
                            # Stop calculation stages                                                           
                            break
                        } else {
                            UpdateGeneralInfo calculation_state "OKcalc" $stageNode
                        }
                    }  

                    # Read CPS_$number_last_step
                    set xp_num_steps {container[@n="Calculation_Data"]/container[@n="CALCULATION_STEP_DATA"]/}
                    append xp_num_steps {value[@n="number_of_calculation_steps"]}                            
                    set num_stepsNode [$stageNode selectNodes $xp_num_steps]                           
                    set num_steps [$num_stepsNode @v]
                    set last_calc_step [expr {$num_steps + 1}]                            
                    if { $last_calc_step < 10 } {
                        set cps_extension "CPS_00$last_calc_step"
                    } else {
                        set cps_extension "CPS_0$last_calc_step"
                    }                                        
                    
                    # Read TOTAL TIME, NUMBER_OF_ACTIVE_ELEMENTS AND OUTPUT MATERIAL POINTS
                    set CPS_A3D_file [file join $A3D_folder "$projectName.$cps_extension"]
                    if { [file exists $CPS_A3D_file] } {
                        set f [open $CPS_A3D_file r]
                        set is_total_time 0
                        set is_number_of_active_elements 0
                        set is_out_num_material_points 0
                        set is_out_material_points 0
                        set out_material_points ""
                        set i 0
                        set fcontents [split [read $f] \n]
                        close $f
                        foreach line $fcontents {
                            if { $line eq {$$TOTAL_TIME} } {
                                set is_total_time 1
                                continue
                            }                            
                            if { $is_total_time } { 
                                set total_time $line
                                set is_total_time 0 
                            }
                            if { $line eq {$$NUMBER_OF_ACTIVE_ELEMENTS} } {
                                set is_number_of_active_elements 1
                                continue
                            }                            
                            if { $is_number_of_active_elements } { 
                                set number_of_active_elements $line 
                                set is_number_of_active_elements 0
                            }
                            if { $line eq {$$OUTPUT_NUMBER_OF_MATERIAL_POINTS} } {
                                set is_out_num_material_points 1
                                continue
                            }                            
                            if { $is_out_num_material_points } { 
                                set out_num_material_points $line
                                set is_out_num_material_points 0 
                            }
                            if { $line eq {$$OUTPUT_MATERIAL_POINTS} } {
                                set is_out_material_points 1
                                continue
                            }                            
                            if { $is_out_material_points } {
                                incr i
                                if { $i <= $out_num_material_points } {
                                    lappend out_material_points $line  
                                    set is_out_material_points 1
                                } else {
                                    set is_out_material_points 0
                                }
                                continue                                 
                            }
                        }                                                                                                              
                    }
					UpdateGeneralInfo cps_extension_after_calc $cps_extension $stageNode
					UpdateGeneralInfo total_time $total_time $stageNode        
                    UpdateGeneralInfo prev_stage_steps $num_steps $stageNode        
                    UpdateGeneralInfo number_of_active_elements $number_of_active_elements $stageNode
					UpdateGeneralInfo out_num_material_points $out_num_material_points $stageNode
					UpdateGeneralInfo out_material_points $out_material_points $stageNode
                    
					
                    set CPS_A3D_file [file join $A3D_folder "$projectName.$cps_extension"]
                    if { [file exists $CPS_A3D_file] } { 
                        file delete -force $CPS_A3D_file
                    }                   
                    set 2dat_file [file join $projectfolder "$projectName-2.dat"]
                    if { [file exists $2dat_file] } { 
                        file delete -force $2dat_file
                    }
					
					}                 
                incr icount_stage    
            } 
            global RunProcInfo
            set RunProcInfo "" 
            set calc_state_List ""           
            set xp {container[@n='stages']/blockdata[@n='stage']}  
            set stagesnames ""
            foreach stageNode [$root selectNodes $xp]  {   
                lappend stagesnames [$stageNode @name]              
                set stagetocalculate 0                
                if { $args eq "local" } {                   
                    set stagetocalculate 1
                } else {                    
                    set stagetocalculate [Anura3D::StageToCalculate [$stageNode @name] $args]
                }                    
                if { $stagetocalculate } {
                    lappend calc_state_List [ReadGeneralInfo calculation_state $stageNode]
                    
                }            
            }
            if { [lsearch $calc_state_List "wrongcalc"] == -1 } { 
                if { $args eq "local" } {
                    if { ![isalive $pid] } {
                        snit_messageBox -parent .gid -message \
                            [= "All the stages has been successfully calculated"]               
                        return
                    }
                } else {
                    snit_messageBox -parent .gid -message \
                        [= "The stages from the 'Stages manager' \nhas been successfully calculated: %s" [join $stagesnames ","]]               
                    return
                }
            } 
        }
    }       
}

proc isalive { pid } {
    return [ gid_cross_platform::process_exists $pid ]
}

proc Anura3D::anura3D_WriteBatchFile { filebatch icount_stage cps_extension gom_extension} {
    GiD_WriteCalculationFile init $filebatch
    GiD_WriteCalculationFile puts "@echo off"
    
    GiD_WriteCalculationFile puts "rem basename = %1"
    GiD_WriteCalculationFile puts "rem modeldirectory = %2"
    GiD_WriteCalculationFile puts "rem problemtypedir = %3"
    GiD_WriteCalculationFile puts "rem A3Dfolder = %4"
    GiD_WriteCalculationFile puts "rem OutputFile = %4\\%1.PROC"
    GiD_WriteCalculationFile puts "rem ErrorFile: %2\\%1.err"
    GiD_WriteCalculationFile puts "if exist %2\\%1.POST.lst del %2\\%1.POST.lst"
    GiD_WriteCalculationFile puts "if exist %4\\%1.BMR del %4\\%1.BMR"
    GiD_WriteCalculationFile puts "if exist %4\\%1.BMS del %4\\%1.BMS"
    GiD_WriteCalculationFile puts "if exist %4\\%1.TST del %4\\%1.TST"
    GiD_WriteCalculationFile puts "if exist %4\\%1_000.RX del %4\\%1_000.RX"        
    GiD_WriteCalculationFile puts "copy %2\\%1-1.dat %2\\%1.$cps_extension"    
    GiD_WriteCalculationFile puts "if exist %4\\%1.$cps_extension del %4\\%1.$cps_extension"
    GiD_WriteCalculationFile puts "copy %2\\%1.$cps_extension %4\\%1.$cps_extension"
    GiD_WriteCalculationFile puts "copy %2\\%1-2.dat %2\\%1.$gom_extension"
    GiD_WriteCalculationFile puts "if exist %4\\%1.$gom_extension del %4\\%1.$gom_extension"    
    GiD_WriteCalculationFile puts "copy %2\\%1.$gom_extension %4\\%1.$gom_extension"      
    GiD_WriteCalculationFile end           
}

proc Anura3D::LoadScripts { } {
    variable problemtype_dir
    # Common scripts
    set script_files [list createGOM.tcl createCPS.tcl \
            stagesmanager.tcl materialparameters.tcl \
            createPSF.tcl createHHBF.tcl watertable.tcl hydraulichead.tcl\
            velocity_in_time.tcl createPVF.tcl createOPD.tcl]
    
    foreach filename $script_files {
        uplevel #0 [list source [file join $problemtype_dir scripts $filename]]
    }
}

### auxiliary procedures to be called from the .bas templates to write some of its parts

proc Anura3D::WriteDamping { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data_by_element([lindex $item 1]) [lindex $item 3]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data_by_element($element_id)] } {
            append result "$data_by_element($element_id) \n"
        } else {
            append result "0 \n"
        }
    }
    return $result
}

proc Anura3D::WriteNumberOfMaterialPointsSP { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data_by_element([lindex $item 1]) [lindex $item 3]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data_by_element($element_id)] } {
            append result "$data_by_element($element_id) 0 \n"
        } else {
            append result "0 0 \n"
        }
    }
    return $result
}

proc Anura3D::WriteNumberOfMaterialPointsDP { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data1_by_element([lindex $item 1]) [lindex $item 3]
                set data2_by_element([lindex $item 1]) [lindex $item 4]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data1_by_element($element_id)] || [info exists data2_by_element($element_id)] } {
            append result "$data1_by_element($element_id) $data2_by_element($element_id) \n"
        } else {
            append result "0 0 \n"
        }
    }
    return $result
}

#to be used by TKWIDGET to adda layer thickness
#e.g.
#QUESTION: your_question
#VALUE: soil_thickness
#TKWIDGET: GidUtils::TkwidgetPickSoilLayer

proc GidUtils::TkwidgetPickSoilLayer { event args } {
    global tkwidgedprivpicknodebuttons
    switch $event {
        INIT {
            lassign $args PARENT current_row_variable GDN STRUCT QUESTION
            upvar $current_row_variable ROW
            set entry ""
            set entry_gridded 0
            foreach item [grid slaves $PARENT -row [expr $ROW-1]] {
                if { [winfo class $item] == "Entry"  || [winfo class $item] == "TEntry" } {
                    #assumed that it is the only entry of this row
                    set entry $item
                    set entry_gridded 1
                    break
                }
            }
            if { $entry == "" } {
                set entry [GidUtils::DarkTrickTryFindUngriddedEntry $PARENT $QUESTION]
            }
            if { $entry != "" } {
                set tkwidgedprivpicknodebuttons($QUESTION) [ttk::button $PARENT.bpicknode$QUESTION \
                        -image [gid_themes::GetImage "dimension_dist.png" small_icons] \
                        -command [list GetThickness $entry]]
                grid $tkwidgedprivpicknodebuttons($QUESTION) -row [expr $ROW-1] -column 2 -sticky w
                grid configure $entry -sticky ew
                if { !$entry_gridded } {
                    grid remove $entry
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
            return ""
        }
        SYNC {
            #lassign $args GDN STRUCT QUESTION
            #DWLocalSetValue $GDN $STRUCT $QUESTION $value
        }
        DEPEND {
            lassign $args GDN STRUCT QUESTION ACTION VALUE
            if { [info exists tkwidgedprivpicknodebuttons($QUESTION)] && \
                     [winfo exists $tkwidgedprivpicknodebuttons($QUESTION)] } {
                if { $ACTION == "HIDE" } {
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                } else {
                    #RESTORE
                    grid $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
        }
        CLOSE {
            array unset tkwidgedprivpicknodebuttons
        }
        default {
            return [list ERROR [_ "Unexpected tkwidget event"]]
        }
    }
    #a tkwidget procedure must return "" if Ok or [list ERROR $description] or [list WARNING $description]
    return ""
}

#Required for the previous function
proc GetThickness {entry} {
    set p1 [GidUtils::GetCoordinates [_ "Enter first point (ESC to leave)"]]
    if { $p1=="" } return ""
    set p2 [GidUtils::GetCoordinates [_ "Enter second point (ESC to leave)"]]
    if { $p2=="" } return ""

    set d [MathUtils::VectorDistance $p1 $p2]

    $entry delete 0 end
    $entry insert end $d

}


#to be used by TKWIDGET to pick soil surface
#e.g.
#QUESTION: your_question
#VALUE: soil_surface (Y coordinate)
#TKWIDGET: GidUtils::TkwidgetSoilElevation

proc GidUtils::TkwidgetSoilElevation { event args } {
    global tkwidgedprivpicknodebuttons
    switch $event {
        INIT {
            lassign $args PARENT current_row_variable GDN STRUCT QUESTION
            upvar $current_row_variable ROW
            set entry ""
            set entry_gridded 0
            foreach item [grid slaves $PARENT -row [expr $ROW-1]] {
                if { [winfo class $item] == "Entry"  || [winfo class $item] == "TEntry" } {
                    #assumed that it is the only entry of this row
                    set entry $item
                    set entry_gridded 1
                    break
                }
            }
            if { $entry == "" } {
                set entry [GidUtils::DarkTrickTryFindUngriddedEntry $PARENT $QUESTION]
            }
            if { $entry != "" } {
                set tkwidgedprivpicknodebuttons($QUESTION) [ttk::button $PARENT.bpicknode$QUESTION \
                        -image [gid_themes::GetImage "point.png" small_icons] \
                        -command [list GetElevation $entry]]
                grid $tkwidgedprivpicknodebuttons($QUESTION) -row [expr $ROW-1] -column 2 -sticky w
                grid configure $entry -sticky ew
                if { !$entry_gridded } {
                    grid remove $entry
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
            return ""
        }
        SYNC {
            #lassign $args GDN STRUCT QUESTION
            #DWLocalSetValue $GDN $STRUCT $QUESTION $value
        }
        DEPEND {
            lassign $args GDN STRUCT QUESTION ACTION VALUE
            if { [info exists tkwidgedprivpicknodebuttons($QUESTION)] && \
                     [winfo exists $tkwidgedprivpicknodebuttons($QUESTION)] } {
                if { $ACTION == "HIDE" } {
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                } else {
                    #RESTORE
                    grid $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
        }
        CLOSE {
            array unset tkwidgedprivpicknodebuttons
        }
        default {
            return [list ERROR [_ "Unexpected tkwidget event"]]
        }
    }
    #a tkwidget procedure must return "" if Ok or [list ERROR $description] or [list WARNING $description]
    return ""
}

#Required for the previous function
proc GetElevation {entry} {
    set p1 [GidUtils::GetCoordinates [_ "Enter first point (ESC to leave)"]]
    if { $p1=="" } return ""
    set Elev [lindex $p1 1]
    $entry delete 0 end
    $entry insert end $Elev
}


######################################################################
#  auxiliary procs invoked from the tree (see .spd xml description)


proc Anura3D::GetMaterialsList { domNode } {
    set xp {//container[@n="materials"]}
    set dom_materials [$domNode selectNodes $xp] 
    set result ""
    foreach dom_matNode $dom_materials {
        foreach dom_material [$dom_matNode childNodes] {
            set name [$dom_material @name]
            if { [lsearch -exact $result $name] == "-1" } {            
                lappend result $name
            }
        }   
    }      
    return [join $result ,]
}

proc Anura3D::EditDatabase { domNode dict boundary_conds } {   
    set has_container ""
    set database materials
    set title [= "User defined"]
    set list_name [$domNode @n]
    
    set xp {ancestor::blockdata[@n="stage"]}
    set stageNode [$domNode selectNodes $xp]
    set stage_name [$stageNode @name]
    
    set xp  [format_xpath {/*/container[@n="stages"]/blockdata[@name=%s]/container[@n=%s]} \
            $stage_name "materials"]
    
    set dom_materials [$domNode selectNodes $xp]     
  
    set primary_level "material"
    if { [dict exists $dict $list_name] } {      
        set xps [format_xpath {/*/container[@n="stages"]/blockdata[@name=%s]/container[@n=%s]/blockdata[@n=%s and @name=%s]} \
                $stage_name "materials" "material" [dict get $dict $list_name]]                
    } else {
        set xps ""
    }
    set domNodes [gid_groups_conds::edit_tree_parts_window -accepted_n "material" \
            -select_only_one 1 $boundary_conds $title $xp $xps]
    set dict ""
    if { [llength $domNodes] } {
        set domNode [lindex $domNodes 0]
        if { [$domNode @n] eq $primary_level } {
            dict set dict $list_name [$domNode @name]
        }
    }
    return [list $dict ""]
}

proc check_dim_points {dim1 dim2 npoints domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set point_path {string(//container[@n="General_data"]/value[@n="NLAYERS"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    set num_points [$domNode selectNodes $point_path]
    if {[string equal $problem_type $dim1]} {
        if {[string equal $num_points $npoints]} {
		return normal}
    }
    if {[string equal $problem_type $dim2]} {
        if {[string equal $num_points $npoints]} {return normal}
    }
    return hidden
}

proc check_dim {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return normal}
    if {[string equal $problem_type $dim2]} {return normal}
    return hidden
}

proc check_points {point domNode} {
    set point_path {string(//container[@n="General_data"]/value[@n="NLAYERS"]/@v)}
    set num_points [$domNode selectNodes $point_path]
    if {[string equal $num_points $point]} {return normal}
    return hidden
}

proc check_dim_surface {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return line}
    if {[string equal $problem_type $dim2]} {return line}
    return line,surface
}

proc check_dim_linesurface_surfacevolume {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return line,surface}
    if {[string equal $problem_type $dim2]} {return line,surface}
    return surface,volume
}

proc check_dim_line_surface {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return line}
    if {[string equal $problem_type $dim2]} {return line}
    return surface
}

proc check_dim_surface_volume {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return surface}
    if {[string equal $problem_type $dim2]} {return surface}
    return volume
}

proc check_dim_prescribed_velocity {dim1 dim2 domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {[string equal $problem_type $dim1]} {return point,line,surface}
    if {[string equal $problem_type $dim2]} {return point,line,surface}
    return point,line,surface,volume
}

proc hide_show_y_n_liq_n {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Liquid"} {return hidden}
    return normal
}

proc hide_show_y_n_dry {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Dry material"} {return hidden}
    return normal
}

proc hide_show_y_n_intr_perm {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Saturated material-fully coupled"} {return normal}
    if {$problem_type == "Unsaturated material-2-phase with suction effect"} {return normal}
    if {$problem_type == "Unsaturated material-3-phase fully coupled"} {return normal}
    return hidden
}

proc hide_show_y_n_liq_y {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Liquid"} {return normal}
    return hidden
}

proc hide_show_y_n_sat_ful {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Saturated material-fully coupled"} {return normal}
    if {$problem_type == "Unsaturated material-2-phase with suction effect"} {return normal}
    if {$problem_type == "Unsaturated material-3-phase fully coupled"} {return normal}
    if {$problem_type == "Liquid"} {return normal}
    return hidden
}

proc hide_show_y_n_un_3 {domNode} {
    set dim_path {string(../value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == "Unsaturated material-3-phase fully coupled"} {return normal}
    return hidden
}

proc hide_show_const_model {flag domNode} {
    set x_path {string(../value[@n="_material_model_solid_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type == "Liquid"} {return hidden}
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc hide_show_const_model_mat_type_y {flag domNode} {
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc hide_show_const_model_mat_type_n {flag domNode} {
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    if {$problem_type == $flag} {return hidden}
    return normal
}

proc hide_show_const_model_z {flag domNode} {
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}
    set dim_type [$domNode selectNodes $dim_path]
    set x_path {string(../value[@n="_material_model_solid_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    if {$dim_type == "2D:plane-strain"} {return hidden}
    if {$dim_type == "2D:Axissymmetric"} {return hidden}
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc hide_show_const_model_multi_undr {flag1 flag2 domNode} {
    set dim_path {string(../value[@n="_material_model_solid_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type == "Liquid"} {return hidden}
    if {$problem_type == $flag1 && $mat_type == "Saturated material-undrained effective stress"} {return normal}
    if {$problem_type == $flag2 && $mat_type == "Saturated material-undrained effective stress"} {return normal}
	return hidden
}

proc hide_show_const_model_multi {flag1 flag2 domNode} {
    set dim_path {string(../value[@n="_material_model_solid_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type == "Liquid"} {return hidden}
    if {$problem_type == $flag1} {return normal}
    if {$problem_type == $flag2} {return normal}
    return hidden
}
proc hide_show_const_model_multi_liq {flag1 flag2 domNode} {
    set dim_path {string(../value[@n="_material_model_liquid_"]/@v)}
    set problem_type [$domNode selectNodes $dim_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type != "Liquid"} {return hidden}
    if {$problem_type == $flag1} {return normal}
    if {$problem_type == $flag2} {return normal}
    return hidden
}

proc hide_show_const_model_liq {flag domNode} {
    set x_path {string(../value[@n="_material_model_liquid_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type != "Liquid"} {return hidden}
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc free_surf_fact {domNode} {
    set x_path {string(../value[@n="detect_free_surface_liquid_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type != "Liquid"} {return hidden}
    if {$problem_type == "Yes"} {return normal}
    return hidden
}

proc is_unsat {domNode} {
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type == "Unsaturated material-2-phase with suction effect"} {return normal}
    if {$mat_type == "Unsaturated material-3-phase fully coupled"} {return normal}
    return hidden
}

proc is_unsat_ret {flag domNode} {
    set x_path {string(../value[@n="_unsat_retention_curve_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type != "Unsaturated material-2-phase with suction effect"} {
        if {$mat_type != "Unsaturated material-3-phase fully coupled"} {return hidden}
    }
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc is_unsat_cond { flag domNode } {
    set x_path {string(../value[@n="_unsat_hydraulic_cond_"]/@v)}
    set problem_type [$domNode selectNodes $x_path]
    set dim_path {string(../../container[@n="_basic"]/value[@n="material_type_"]/@v)}
    set mat_type [$domNode selectNodes $dim_path]
    if {$mat_type != "Unsaturated material-2-phase with suction effect"} {
        if {$mat_type != "Unsaturated material-3-phase fully coupled"} {return hidden}
    }
    if {$problem_type == $flag} {return normal}
    return hidden
}

proc find_material_id { material_name stageNode } {
    set xp {container[@n="materials"]/blockdata}    
    set int 1
    foreach gNode [$stageNode selectNodes $xp] {                       
        if {[$gNode @name] eq $material_name} {
            return $int
        }            
        incr int
    }
    # If no matching blockdata element is found, return an empty string.
    return ""
}

proc activate_initial_stage { domNode dict bc } {
    set stagename [stage_name $domNode]
    set stageList [stage_list]
    set initial_stage 0
    set state hidden
    if { $stagename eq [lindex $stageList 0]} {
        set initial_stage 1 
        set state normal 
    }
    return $state     
}

proc activate_not_initial_stage { domNode } {
    set stagename [stage_name $domNode]
    set stageList [stage_list]
    set not_initial_stage 0
    set state hidden
    if { $stagename != [lindex $stageList 0]} {
        set not_initial_stage 1 
        set state normal 
    }
    return $state     
}

proc stage_name { domNode } {         
    set xp {ancestor::blockdata[@n="stage"]}
    set stageNode [$domNode selectNodes $xp]  
    return [$stageNode @name]
}

proc stage_list {} {
    set root [$gid_groups_conds::doc documentElement]
    set stageList ""
    set xp {container[@n='stages']/blockdata[@n='stage']}      
    foreach stageNode [$root selectNodes $xp]  { 
        lappend stageList [$stageNode @name]  
    }
    return $stageList 
}

proc Anura3D::Bitmaps { dir { type "DEFAULT INSIDELEFT"} } {
    global BitmapsNames BitmapsCommands BitmapsHelp
    global ProgramName VersionNumber Priv
    
    set BitmapsNames(0) [list images/file_tree.png images/post_toolbar.png \
        --- images/stage_22.png images/configure_22.png \
        --- images/mesh.png images/calculate_22.png]
    
    set BitmapsCommands(0) [list \
            {-np- gid_groups_conds::open_conditions menu} \
            {-np- gid_groups_conds::open_conditions menu -select_xpath {/*/container[@n="General_data"]}} \
            {} {-np- gid_groups_conds::open_conditions menu_or_any -select_xpath {/*/container[@n="stages"]}} \
            {-np- Anura3D::StagesManager} {} {-np-  GiD_Process MEscape Meshing generate} {-np- Anura3D::Calculate local}]   

    set BitmapsHelp(0) [list [= "Data tree view"] [= "Define postprocess visualization software"] \
            "" [= "Assign stages"] [= "Stages manager"] "" [= "Generate mesh"] \
        [= "Start calculation process"]  [= "Cancel calculation process"]] 
    

    # prefix values:
    #          Pre        Only active in the preprocessor
    #          Post       Only active in the postprocessor
    #          PrePost    Active Always

    set prefix Pre
    set Priv(toolbarwin) [CreateOtherBitmaps Anura3DBar "Anura3D toolbar" \
            BitmapsNames BitmapsCommands \
            BitmapsHelp $dir [list Anura3D::Bitmaps $dir] $type $prefix]
    AddNewToolbar "$ProgramName bar" ${prefix}BarWindowGeom \
        [list Anura3D::Bitmaps $dir] [= "%s bar" $ProgramName]
}

proc Anura3D::ReadGeneralInfo { type domNode } {      
    set dict [split [$domNode @generalinfo] ,]       
    return [dict get $dict $type]    
}

proc Anura3D::UpdateGeneralInfo { type value domNode } {
    set dict [split [$domNode @generalinfo] ,]        
    dict set dict $type $value   
    set dict [join $dict ","]    
    
    gid_groups_conds::setAttributes [gid_groups_conds::nice_xpath $domNode] \
        [list generalinfo $dict]           
}


proc apply_hydraulic_head { domNode } {    
    set root [$gid_groups_conds::doc documentElement]
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}    
    set dim_type [$root selectNodes $dim_path]
    
    set apply 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/value[@n="Minimum_coordinates_hydraulic_head_2D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin 
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/value[@n="Maximum_coordinates_hydraulic_head_2D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax   
        if { $xmin != $xmax && $ymin != $ymax } {
            set apply 1
        }     
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/value[@n="Minimum_coordinates_hydraulic_head_3D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin zmin
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Hydraulic_head"]/value[@n="Maximum_coordinates_hydraulic_head_3D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax zmax
        if { $xmin != $xmax && $ymin != $ymax && $zmin != $zmax } {
            set apply 1
        }
    }
    return $apply
}

proc apply_seepage_face { domNode } {
    set root [$gid_groups_conds::doc documentElement]
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}    
    set dim_type [$root selectNodes $dim_path]
    
    set apply 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/value[@n="Minimum_coordinates_seepage_face_2D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin 
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/value[@n="Maximum_coordinates_seepage_face_2D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax   
        if { $xmin != $xmax && $ymin != $ymax } {
            set apply 1
        }     
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/value[@n="Minimum_coordinates_seepage_face_3D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin zmin
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Seepage_face"]/value[@n="Maximum_coordinates_seepage_face_3D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax zmax
        if { $xmin != $xmax && $ymin != $ymax && $zmin != $zmax } {
            set apply 1
        }
    }
    return $apply
}

proc apply_infiltration { domNode } {
    set root [$gid_groups_conds::doc documentElement]
    set dim_path {string(//container[@n="General_data"]/value[@n="NDIM"]/@v)}    
    set dim_type [$root selectNodes $dim_path]
    
    set apply 0
    if {$dim_type == "2D:plane-strain" || $dim_type == "2D:Axissymmetric"} {                
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/value[@n="Minimum_coordinates_infiltration_2D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin 
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/value[@n="Maximum_coordinates_infiltration_2D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax   
        if { $xmin != $xmax && $ymin != $ymax } {
            set apply 1
        }     
    } elseif {$dim_type == "3D" || $dim_type == "3D:Axissymmetric"} {
        set min_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/value[@n="Minimum_coordinates_infiltration_3D"]/@v)}]
        set min_coords [split $min_coords ","]
        lassign $min_coords xmin ymin zmin
        set max_coords [$domNode selectNodes {string(../../container[@n="BC"]/container[@n="Hydraulic_Conditions"]/container[@n="Infiltration"]/value[@n="Maximum_coordinates_infiltration_3D"]/@v)}]
        set max_coords [split $max_coords ","]
        lassign $max_coords xmax ymax zmax
        if { $xmin != $xmax && $ymin != $ymax && $zmin != $zmax } {
            set apply 1
        }
    }
    return $apply
}

proc apply_contact { domNode } {  
    set xp {../../container[@n="Contact"]/condition[@n="Contact_properties"]/group}
    set apply 0
    if { [$domNode selectNodes $xp] != "" } {
        set apply 1
    }
    
    return $apply
}

proc apply_initial_velo { domNode } {
    set xp {../../container[@n="Initial_cond"]/condition[@n="Initial_MP_velocity"]/group}    
    set apply 0
    if { [$domNode selectNodes $xp] != "" } { 
        set apply 1 
    }   
    
    return $apply
}

proc post_visualization_software { $post_view_softwareNode } {
    set xp {container[@n="General_data"]/container[@n="POSTPROCESS_VISUALIZATION_SOFTWARE"]}
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set post_view_softwareNode [$root selectNodes $xp]
    set paraview [$post_view_softwareNode selectNodes {string(value[@n="Paraview"]/@v)}]
    set GiD_ASCII [$post_view_softwareNode selectNodes {string(value[@n="GiD_ASCII"]/@v)}]
    set GiD_Binary [$post_view_softwareNode selectNodes {string(value[@n="GiD_Binary"]/@v)}]

    if { $paraview && !$GiD_ASCII && !$GiD_Binary } {
        return "ParaView Output"       
    } 
    if { $GiD_ASCII } {
        if { !$paraview } {
            return "GiD Output"
        } else {
            return "ParaView and GiD Output"
        }         
    } 
    if { $GiD_Binary } {
        if { !$paraview } {
            return "GiD Output"
        } else {
            return "ParaView and GiD Output"
        }                 
    }
}


