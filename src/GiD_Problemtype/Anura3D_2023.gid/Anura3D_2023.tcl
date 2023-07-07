proc InitGIDProject { dir } {

    global _dir
    set _dir $dir
    set anura_version "Anura3D v2023"
    set GiDVersionRequired "14.0"

    Anura3D::SetDir $dir ;#store to use it later
    Anura3D::LoadScripts
    GidUtils::OpenWindow CUSTOMLIB

	# disclaimer at start-up
	set self_close 0
    GidUtils::Splash [file join $_dir images about_anura3d.png] .splash $self_close [list $anura_version 445 10]

	# check minimal GiD version
	set message_required "This $anura_version interface is developed for GiD $GiDVersionRequired or later. \n \n It is advised to update your GiD software."
	set title_required "Anura3D - Required GiD Version"
    if { [GidUtils::VersionCmp $GiDVersionRequired] < 0 } { tk_messageBox -title $title_required -message $message_required -icon warning -type ok }

	# add ANURA3D menu
	GiDMenu::Create "Anura3D" "PRE" 5 =
    #GiDMenu::InsertOption "Anura3D" [list "Materials"] 0 PRE "GidOpenMaterials" "Control-A" "[file join $dir images icon_material.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Material Point Specification"] 1 PRE "GidOpenConditions Material_Point_Specification" "Control-B" "[file join $dir images icon_mps.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Fixities"] 2 PRE "GidOpenConditions Fixities" "Control-C" "[file join $dir images icon_fixity.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Remove Fixities"] 3 PRE "GidOpenConditions Remove_Fixities" "Control-D" "[file join $dir images icon_remove.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Loading Conditions"] 4 PRE "GidOpenConditions Loading_Conditions" "Control-E" "[file join $dir images icon_load.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Soil Surface"] 5 PRE "GidOpenConditions Soil_Surface" "Control-F" "[file join $dir images icon_soil_surface.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Phreatic Surface"] 6 PRE "GidOpenConditions Phreatic_Surface" "Control-G" "[file join $dir images icon_phreatic.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Prescribed Velocities"] 7 PRE "GidOpenConditions Prescribed_Velocities" "Control-H" "[file join $dir images icon_velocity.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Initial Conditions"] 8 PRE "GidOpenConditions Initial_Conditions" "Control-I" "[file join $dir images icon_initial.png]" replace =
	#GiDMenu::InsertOption "Anura3D" [list "Damping Conditions"] 9 PRE "GidOpenConditions Damping_Conditions" "Control-J" "[file join $dir images icon_damp.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Contact Properties"] 10 PRE "GidOpenConditions Contact_Properties" "Control-L" "[file join $dir images icon_contact.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Excavation"] 11 PRE "GidOpenConditions Excavation" "Control-M" "[file join $dir images icon_excavation.png]" replace =
    ## GiDMenu::InsertOption "Anura3D" [list "Construction"] 10 PRE "GidOpenConditions Construction" "Control-K" "[file join $dir images icon_construction.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Reaction Forces"] 12 PRE "GidOpenConditions Reaction_Forces" "Control-N" "[file join $dir images icon_reactionforce.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Absorbing Boundaries"] 13 PRE "GidOpenConditions Absorbing_Boundaries" "Control-O" "[file join $dir images icon_absorbing.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Moving Mesh"] 14 PRE "GidOpenConditions Moving_Mesh" "Control-P" "[file join $dir images icon_movingmesh.png]" replace =
	#GiDMenu::InsertOption "Anura3D" [list "Hydraulic Boundary Conditions"] 15 PRE "GidOpenConditions Hydraulic_Boundary_Conditions" "Control-Q" "[file join $dir images icon_hydraulicBC.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "---"] 88 PRE "" "" "" replace =
    ## GiDMenu::InsertOption "Anura3D" [list "Output Data"] 89 PRE "GidOpenProblemData Select_Output_Data" "Control-X" "[file join $dir images icon_output.png]" replace =
    ## GiDMenu::InsertOption "Anura3D" [list "---"] 90 PRE "" "" "" replace =
    #GiDMenu::InsertOption "Anura3D" [list "Calculation Data"] 0 PRE "GidOpenProblemData Calculation_Data" "Control-A" "[file join $dir images cps.png]" replace =
    #GiDMenu::InsertOption "Anura3D" [list "---"] 1 PRE "" "" "" replace =
    GiDMenu::InsertOption "Anura3D" [list "Generate Anura3D Files"] 0 PRE "Anura3D::Calculate" "Control-B" "[file join $dir images generate.png]" replace =
    GiDMenu::InsertOption "Anura3D" [list "---"] 1 PRE "" "" "" replace =
	GiDMenu::InsertOption "Anura3D" [list "Tutorial Manual..."] 2 PRE "Anura3D::Tutorial" "" "[file join $dir images tutorial.png]" replace =
	GiDMenu::InsertOption "Anura3D" [list "Scientific Manual..."] 3 PRE "Anura3D::Scientific" "" "[file join $dir images scientific.png]" replace =
	GiDMenu::InsertOption "Anura3D" [list "Verification Manual..."] 4 PRE "Anura3D::Verification" "" "[file join $dir images verification.png]" replace =
    GiDMenu::InsertOption "Anura3D" [list "Disclaimer..."] 5 PRE "Anura3D::Disclaimer" "" "[file join $dir images icon_anura3d.png]" replace =
	GiDMenu::InsertOption "Anura3D" [list "About..."] 6 PRE "Anura3D::About" "" "[file join $dir images icon_anura3d.png]" replace =

	# simplify and adaptd HELP menu
	#GiDMenu::RemoveOption "Help" [list "Customization help"] "PRE" _
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
    GiDMenu::UpdateMenus

    # open window tree
    GidUtils::OpenWindow CUSTOMLIB

}

proc EndGIDProject { } {

}

namespace eval Anura3D {
}

namespace eval ANURA3D_2023 {
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
    set TestDoc [file join $_dir doc "TutorialManual_2023.pdf"]
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

proc Anura3D::About { } {

    tk_messageBox -title "Anura3D - About" -message "Anura3D Version 2023\n\nFor more information about Anura3D, check the website\n\nhttp://www.anura3D.com" -icon info

}

proc Anura3D::Disclaimer { } {

set answer [tk_messageBox -title "Anura3D - Disclaimer" -message "Copyright (C) 2020  Members of the Anura3D MPM Research Community

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>

Anura3D MPM Research Community
E-mail:	info@anura3D.com
Web: 	www.anura3D.com" -icon info -type ok]

}

#proc Anura3D::CreateBatch { channel } {

#    set project_path [GiD_Info Project ModelName]
#    set model_name [file tail $project_path]
#    set exe_name [GiD_Info Project ProblemType]
#    GiD_File fprintf -nonewline $channel "\""
#    GiD_File fprintf -nonewline $channel $project_path
#    GiD_File fprintf -nonewline $channel ".A3D\\"
#    GiD_File fprintf -nonewline $channel $exe_name
#    GiD_File fprintf -nonewline $channel ".exe\" "
#    GiD_File fprintf -nonewline $channel "\""
#    GiD_File fprintf -nonewline $channel $project_path
#    GiD_File fprintf -nonewline $channel ".A3D\\"
#    GiD_File fprintf -nonewline $channel $model_name
#    GiD_File fprintf $channel "\""
#    GiD_File fprintf $channel "PAUSE"

#}

proc Anura3D::Calculate { } {
	set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
	set calculation_name $project_path
	set CPS_name $project_path
	set GOM_name $project_path
	set OPD_name $project_path
	append calculation_name ".gid" "/" $model_name ".dat"
	append CPS_name ".gid" "/" $model_name "-1.dat"
	append GOM_name ".gid" "/" $model_name "-2.dat"
	append OPD_name ".gid" "/" $model_name "-3.dat"
	
    set answer [tk_messageBox -title "Anura3D - Generate Anura3D Files" -message "Before creating the Anura3D input files, make sure that: \n \n \
 - Materials are assigned \n  - Boundary conditions are defined \n  - Material points are specified \n  - Calculation parameters are defined \n \
 - Mesh is (re-)generated.\n\nSave the project files and make sure that the project name\ndoes NOT contain any spaces!" -icon warning -type okcancel]

    switch -- $answer {
      cancel return
      ok {
	  GiD_Process Mescape Utilities Calculate
	  Anura3D::WriteCalculationFile $calculation_name
      Anura3D::WriteCalculationFile_CPS $CPS_name
	  Anura3D::WriteCalculationFile_GOM $GOM_name
	  Anura3D::WriteCalculationFile_OPD $OPD_name
	  }
    }
}


proc Anura3D::LoadScripts { } {
    variable problemtype_dir
    # Common scripts
    set script_files [list calculate.tcl createGOM.tcl createCPS.tcl createOPD.tcl]
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
proc ANURA3D_2023::GetMaterialsList { domNode } {
    set x_path {//container[@n="materials"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set result [list]
    foreach dom_material [$dom_materials childNodes] {
        set name [$dom_material @name]
        lappend result $name
    }
    return [join $result ,]
}

proc ANURA3D_2023::EditDatabaseListDirect {domNode dict boundary_conds } {
    set has_container ""
    set database materials
    set title [= "User defined"]
    set list_name [$domNode @n]
    set x_path {//container[@n="materials"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set primary_level material
    if { [dict exists $dict $list_name] } {
        set xps $x_path
        append xps [format_xpath {/blockdata[@n=%s and @name=%s]} $primary_level [dict get $dict $list_name]]
    } else {
        set xps ""
    }
    set domNodes [gid_groups_conds::edit_tree_parts_window -accepted_n $primary_level -select_only_one 1 $boundary_conds $title $x_path $xps]
    set dict ""
    if { [llength $domNodes] } {
        set domNode [lindex $domNodes 0]
        if { [$domNode @n] == $primary_level } {
            dict set dict $list_name [$domNode @name]
        }
    }
    return [list $dict ""]
}

proc check_dim_points {dim1 dim2 npoints domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set point_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
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
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  if {[string equal $problem_type $dim1]} {return normal}
	  if {[string equal $problem_type $dim2]} {return normal}
      return hidden
}

proc check_points {point domNode} {
      set point_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
	  set num_points [$domNode selectNodes $point_path]
	  if {[string equal $num_points $point]} {return normal}
      return hidden
}

proc check_dim_surface {dim1 dim2 domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  if {[string equal $problem_type $dim1]} {return line}
	  if {[string equal $problem_type $dim2]} {return line}
      return line,surface
}

proc check_dim_linesurface_surfacevolume {dim1 dim2 domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  if {[string equal $problem_type $dim1]} {return line,surface}
	  if {[string equal $problem_type $dim2]} {return line,surface}
      return surface,volume
}

proc check_dim_line_surface {dim1 dim2 domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  if {[string equal $problem_type $dim1]} {return line}
	  if {[string equal $problem_type $dim2]} {return line}
      return surface
}

proc check_dim_surface_volume {dim1 dim2 domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  if {[string equal $problem_type $dim1]} {return surface}
	  if {[string equal $problem_type $dim2]} {return surface}
      return volume
}

proc check_dim_prescribed_velocity {dim1 dim2 domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
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
	  set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set dim_type [$domNode selectNodes $dim_path]
      set x_path {string(../value[@n="_material_model_solid_"]/@v)}
	  set problem_type [$domNode selectNodes $x_path]
	  if {$dim_type == "2D:plane-strain"} {return hidden}
	  if {$dim_type == "2D:Axissymmetric"} {return hidden}
	  if {$problem_type == $flag} {return normal}
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

proc is_unsat_cond {flag domNode} {
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

proc find_material_id {material_name root} {
    set xp [format_xpath {container[@n="materials"]/blockdata}]
    set list [$root selectNodes $xp]
    set int 1
    foreach gNode $list {	
	    set type [$gNode getAttribute "name"]	
	    if {$type eq $material_name} {
	      return $int
	    }	    
	    set int [expr $int + 1]    
    }
    # If no matching blockdata element is found, return an empty string.
    return ""
}
