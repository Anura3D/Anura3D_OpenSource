#######################################################
# WATER TABLE FOR PHREATIC SURFACE
#######################################################

package require fulltktree 
package require base64

namespace eval WaterTable {                   
}

proc WaterTable::create_table { domNode } {
    global Priv     
    
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]           
    
    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    }          
    
    set units_mesh [gid_groups_conds::give_mesh_unit]
    
    set w .gid.ttktable
    destroy $w
    
    set xyList ""    
    set containerNode [$domNode selectNodes {container}]
    if { $containerNode != "" } {               
        foreach iNode [$containerNode selectNodes {value}] {
            lappend xyList [$iNode @v]
        }
    } 
    
    set w [dialogwin_snit $w -title [= "Add water table"] \
            -okname [_ "Ok"] -style ridgeframe -grab 1]
    set f [$w giveframe]
     
    set f1 [ttk::panedwindow $f.f1 -orient horizontal]
    
    set columns [list \
            [list 10 "x \[$units_mesh\]" left text 0 ] \
            [list 8 "y \[$units_mesh\]" left text 0] \
            ]
    
    set f11 [ttk::labelframe $f1.f1 -text [_ "Coordinates entry"]]
    set fgraph [ttk::labelframe $f11.fgraph -text [_ "Graph"]]
    
    $w set_uservar_value fgraph $fgraph 
    
    set tree [fulltktree $f1.left -columns $columns -height 420 -width 100 \
            -expand 1 -bd 1 \
            -selectmode extended -expand 1 -bd 1 -relief solid \
            -showlines 0 -sensitive_cols all \
            -contextualhandler [list WaterTable::contextualmenu $w] \
            -deletehandler [list WaterTable::delete $w $f1.left selection]]
        
    $tree notify bind [$tree givetreectrl] <Header-invoke> ""
    $w set_uservar_value tree $tree  
    
    ttk::label $f11.label_x -text "x \[$units_mesh\]"
    ttk::entry $f11.entry_x -textvariable \
        [$w give_uservar entry_x] -width 22
    $w set_uservar_value entry_x "0.0"
    
    ttk::label $f11.label_y -text "y \[$units_mesh\]"
    ttk::entry $f11.entry_y -textvariable \
        [$w give_uservar entry_y] -width 22
    $w set_uservar_value entry_y "0.0"
    
    ttk::button $f11.b1 -text [_ "Add at end"] -command \
        [list WaterTable::add $w $tree end $f11.data] -width 20
    ttk::button $f11.b2 -text [_ "Add before"] -command \
        [list WaterTable::add $w $tree prev $f11.data] -width 20
    ttk::button $f11.b3 -text [_ "Clear all"] -command \
        [list WaterTable::delete $w $tree all] -width 20
    
    lassign [list 0 0] row col
    
    $w set_uservar_value interpolator_unit $units_mesh    
            
    $f1 add $tree 
    $f1 add $f11
       
    grid $f11.label_x -row 0 -column 0 -sticky w -padx 3 -pady 3
    grid $f11.entry_x -row 0 -column 1 -sticky w -padx 3 -pady 3
    grid $f11.label_y -row 1 -column 0 -sticky w -padx 3 -pady 3
    grid $f11.entry_y -row 1 -column 1 -sticky w -padx 3 -pady 3
        
    grid $f11.b1 $f11.b2 $f11.b3 -sticky news -padx 3 -pady 3
    
    grid $f1 -sticky news -padx 2 -pady 2 
          
    package require customLib_utils::img
   
    package require gid_graph
    
    set points [list]
    set y_text [list]
    foreach col [lrange $columns 1 end] {
        lappend y_text [lindex $col 1]
    }
    set plot_command [list GidGraph::InFrame2 -points $points \
            -title  [_ "Function graph"] \
            -x_text "x-coord" -y_text "y-coord"]
    lappend plot_command $fgraph
    $w set_uservar_value draw_graphs $plot_command
    
    grid $fgraph -columnspan 4 -sticky news
    
    grid columnconfigure $fgraph 0 -weight 1
    grid rowconfigure $fgraph 0 -weight 1
    
    grid columnconfigure $f 0 -weight 1
    grid rowconfigure $f 0 -weight 1
    
    foreach xy $xyList {
        $tree insert end $xy
    }
    
    WaterTable::actualize_graph $w $tree
    
    bind $w <Return> [list $w invokecancel]
    set action [$w createwindow]    
    
    while 1 {
        switch -- $action {             
            -1 - 0 {
                # Cancel                                              
                catch { destroy $w }   
                return ""
            }
            1 {
                # Ok    
                if { ![llength [$tree item children 0]] } {
                    snit_messageBox -parent $w -message \
                        [= "There are no points to define the water profile"]
                } elseif { [llength [$tree item children 0]] < 2 } {
                    snit_messageBox -parent $w -message \
                        [= "There are not enough points to define the water profile"]
                }
                WaterTable::save_data $w $tree $domNode 
                catch { destroy $w }
                return
            }
        }        
    }
}

proc WaterTable::add { w tree where f } {

    set integer_fields ""
    set max_col_by_line 3
    
    set fields [list entry_x entry_y]
    
    set max_col_by_line 5           
    
    lassign [list 0 0] row col
    foreach i $fields {
        set v [$w give_uservar_value $i]
        if { ![string is double -strict $v] } {
            set txt $v            
            set t [_ "Error in field %s. It must be a number" $txt]
            snit_messageBox -parent $w -message $t
            return
        }
        if { $col > $max_col_by_line } {
            set col 0
            incr row
        }
        incr col 2
    }

    set pnts ""
    foreach i $fields { 
        set $i [$w give_uservar_value $i]
    }
   
    set pnt [list]
    foreach i $fields { 
        lappend pnt [set $i]
    }            
    lappend pnts $pnt        
           
    switch $where {
        prev {
            set sibling [lindex [$tree selection get] 0]
            if { $sibling eq "" } {
                snit_messageBox -parent $w -message \
                    [_ "It is necessary to select one row in the table"]
                return
            }
        }
        end {
            set sibling 0
        }
    }
    if { $where eq "end" } {
        set pnt0 [$tree item text last]
        if { $pnt0 eq [lindex $pnts 0] } {
            set ret [snit_messageBox -type okcancel -default ok \
                    -message [_ "Row is equal to the previous one. Are you sure to insert?"] \
                    -parent $w]
            if { $ret == "cancel" } { return }
        }
    }
    foreach i $pnts {
        $tree insert $where $i $sibling
    }
    set itemList [$tree range 0 end]
    if { [llength [$tree range 0 end]] > 2 } {
        WaterTable::actualize_graph $w $tree
    }
}

proc WaterTable::actualize_graph { w tree } {
    set tree_item_children [$tree item children 0]
    set num_rows [llength $tree_item_children]   
    if { $num_rows } {       
        set points [list]
        foreach item $tree_item_children {
            lappend points [$tree item text $item]                
        }              
        set cmd_graph [$w give_uservar_value draw_graphs]
        set pos_points [lsearch $cmd_graph -points]
        if { $pos_points != -1 } {
            incr pos_points
            lset cmd_graph $pos_points $points
        }
        set fgraph [$w give_uservar_value fgraph]
        grid $fgraph -columnspan 4 -sticky news
    
        grid columnconfigure $fgraph 0 -weight 1
        grid rowconfigure $fgraph 0 -weight 1
        eval $cmd_graph             
    } 
} 

proc WaterTable::contextualmenu { w tree items x y } {

    catch { destroy $tree._cmenu }
    set menu [menu $tree._cmenu -tearoff 0]

    $menu add command -label [_ "Select all"] \
        -command [list WaterTable::cutcopypaste $w $tree selectall]
    $menu add separator
    $menu add command -label [_ "Copy"] \
        -command [list WaterTable::cutcopypaste $w $tree copy]
    $menu add command -label [_ "Paste"] \
        -command [list WaterTable::cutcopypaste $w $tree paste]
    $menu add command -label [_ "Paste end"] \
        -command [list WaterTable::cutcopypaste $w $tree paste_end]
    $menu add separator

    $menu add command -label [_ "Delete"] \
        -command [list WaterTable::delete $w $tree selection]
    
    $menu add command -label [_ "Delete all"] \
        -command [list WaterTable::delete $w $tree all]
    
    tk_popup $menu $x $y
}

proc WaterTable::cutcopypaste { w tree what } {

    switch $what {
        selectall {
            $tree selection add all
        }
        cut {
            set data ""
            foreach i [$tree selection get] {
                lappend data [join [$tree item text $i] \t]
                $tree item delete $i
            }
            clipboard clear
            clipboard append [join $data \n]
        }
        copy {
            set tree_selection [lsort -integer [$tree selection get]]
            set data ""
            foreach i $tree_selection {
                lappend data [join [$tree item text $i] \t]
            }
            clipboard clear
            clipboard append [join $data \n]
        }
        paste {
            if { [$tree item count] <= 1 } {
                #if table is empty paste at the end
                return [WaterTable::cutcopypaste $w $tree paste_end]
            }
            set err [catch { clipboard get } data]
            if { $err } { 
                return
            }
            set idx [lindex [$tree selection get] 0]
            if { $idx eq "" } {
                return
            }            
            set num_columns [llength [$tree cget -columns]]
            foreach i [split $data \n] {
                set l [lrange [split $i \t] 0 end]
                if { [llength $l] != $num_columns } {
                    continue
                }
                $tree insert prev $l $idx
            }
            WaterTable::actualize_graph $w $tree
        }
        paste_end {
            set err [catch { clipboard get } data]
            if { $err } { 
                return
            }
            set num_columns [llength [$tree cget -columns]]
            foreach i [split $data \n] {
                set l [lrange [split $i \t] 0 end]
                if { [llength $l] != $num_columns } {
                    continue
                }
                $tree insert end $l
            }
            WaterTable::actualize_graph $w $tree
        }
    }

}

proc WaterTable::value_type { w f button } {
    eval destroy [winfo children $f]

    if { [string index $pn 0] ne "\u0394" } {
        set delta_pn "\u0394$pn"
    } else {
        set delta_pn "\u0394($pn)"
    }

    set max_col_by_line 3
    set fnames [cu::commalist_to_list $fnames]
    set fields [list x]
    set names [list [list $pn:]]
    set i_graph 1
    foreach fname $fnames {
        lappend fields y$i_graph
        lappend names $fname:
        incr i_graph
    }
    set max_col_by_line 5       
    
    lassign [list 0 0] row col
    foreach i $fields j $names {
        set l [cu::nicelabel $f.l${row}_${col} -width 9]
        foreach {txt tag} $j {
            $l insert end $txt $tag
        }       
        set default 0.0 
        
        set e [ttk::entry $f.e${row}_${col} -textvariable [$w give_uservar $i $default] -width 8]
        bind $e <Return> "$button invoke; tk::TabToWindow $f.e0_0; break"
        if { $col > $max_col_by_line } {
            set col 0
            incr row
        }
        grid $l -row $row -column $col -sticky w -padx 2 -pady 2
        grid $e -row $row -column [expr {$col+1}] -sticky w -padx 2 -pady 2
        incr col 2
    }

    grid columnconfigure $f 3 -weight 1

    tk::TabToWindow $f.e0_0
}

proc WaterTable::delete { w tree what args } {

    switch $what {
        all {
            set ret [snit_messageBox -type okcancel -default ok \
                    -message [_ "Are you sure to delete all points?"] \
                    -parent $w]
            if { $ret == "cancel" } { 
                return 
            } else {
                $tree item delete all   
                set fgraph [$w give_uservar_value fgraph]
                grid remove $fgraph
            }
        }
        selection {
            foreach i [$tree selection get] {
                $tree item delete $i
            }    
            WaterTable::actualize_graph $w $tree     
        }               
    }   
}

proc WaterTable::save_data { w tree domNode } {
         
    foreach node [list blockdata container value edit_command dependencies] {
        dom createNodeCmd element $node
    }
    
    set containerNode [$domNode selectNodes {container}]
    if { $containerNode == "" } {
        $domNode appendChildTag container \
            [list attributes() n "water_table_data" \
                pn "Water table data" state hidden]
        set containerNode [$domNode selectNodes {container}]
    } else {
        foreach iNode [$containerNode selectNodes {value}] {
            $iNode delete
        }
    }
    
    if {[llength [$tree item range 0 end]] > 1 } {
        foreach id [$tree item range 0 end] {
            set x [$tree item text $id 0]
            set y [$tree item text $id 1]
            if {[string is double -strict $x] && [string is double -strict $y]} {
                $containerNode appendChildTag value \
                    [list attributes() n "x y" pn "x y" v "$x $y" state hidden]
            }
        }  
    }  
}