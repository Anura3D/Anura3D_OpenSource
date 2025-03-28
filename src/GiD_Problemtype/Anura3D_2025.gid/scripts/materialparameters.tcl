#######################################################
# MATERIAL PARAMETERS AND STATE VARIABLES WINDOW 
#######################################################

namespace eval MaterialParameters {                   
 
}

proc MaterialParameters::mat_params_window { domNode } {  
    global Priv       
    
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]           
    
    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    }          
    
    set xp {string(value[@n="number_mat_params"]/@v)}
    set number_mat_params [$domNode selectNodes $xp]   
    
    set xp {string(value[@n="number_state_var"]/@v)}
    set number_state_var [$domNode selectNodes $xp]
    
    set items [list $number_mat_params $number_state_var]
    
    set max [tcl::mathfunc::max {*}$items]
    
    set height [expr {$max*30}]
    
    set w .gid.ttknote
    destroy $w
    
    set w [dialogwin_snit $w -title [= "Define material parameters and state variables"] \
            -okname [_ "Ok"] -style ridgeframe -grab 1 -geometry "350x350"]
    set f [$w giveframe]
    
    ttk::frame $f.framemat  
    pack $f.framemat -expand yes -fill both -side top    
    
    ttk::notebook $f.framemat.note
    pack $f.framemat.note -fill both -expand 1 -padx 2 -pady 3
    ttk::notebook::enableTraversal $f.framemat.note
    
    # First pane    
    set frmatparam [ttk::frame $f.framemat.note.frmatparam \
            -borderwidth 1 -relief solid] 
    pack $frmatparam -expand yes -fill both -side top
    
    canvas $frmatparam.c -width 350 -height 350 -xscrollcommand "$frmatparam.xscroll set" \
        -yscrollcommand "$frmatparam.yscroll set"
    ttk::scrollbar $frmatparam.xscroll -orient horizontal -command "$frmatparam.c xview"
    ttk::scrollbar $frmatparam.yscroll -command "$frmatparam.c yview"
    pack $frmatparam.xscroll -side bottom -fill x
    pack $frmatparam.yscroll -side right -fill y
    pack $frmatparam.c -expand yes -fill both -side top
    
    # create frame with widgets
    set frmat [ttk::frame $frmatparam.c.frWidgets \
            -borderwidth 1 -relief solid -width 340 -height $height]
    
    for {set i 1} {$i <= $number_mat_params} {incr i} {
        set mat$i [ttk::label $frmat.lb$i \
                -text [= "Material parameter solid $i"]]
        set matval$i [ttk::entry $frmat.en$i \
                -textvariable "[$w give_uservar matval$i]"]      
        grid $frmat.lb$i -padx 2 -pady 2 -row $i -column 0
        grid $frmat.en$i -padx 2 -pady 2 -row $i -column 1     
    }
    
    $frmatparam.c create window 0 0 -anchor nw -window $frmat
    $frmatparam.c configure -scrollregion [$frmatparam.c bbox all]
    
    $f.framemat.note add $frmatparam -text [= "Material parameters solid"] -underline 0 -padding 2
          
    # Second pane
    set frstatevar [ttk::frame $f.framemat.note.frstatevar \
            -borderwidth 1 -relief solid] 
    pack $frstatevar -expand yes -fill both -side top
    
    canvas $frstatevar.c -width 350 -height 350 -xscrollcommand "$frstatevar.xscroll set" \
        -yscrollcommand "$frstatevar.yscroll set"
    ttk::scrollbar $frstatevar.xscroll -orient horizontal -command "$frstatevar.c xview"
    ttk::scrollbar $frstatevar.yscroll -command "$frstatevar.c yview"
    pack $frstatevar.xscroll -side bottom -fill x
    pack $frstatevar.yscroll -side right -fill y
    pack $frstatevar.c -expand yes -fill both -side top
    
    # create frame with widgets
    set frstate [ttk::frame $frstatevar.c.frWidgets \
            -borderwidth 1 -relief solid -width 340 -height $height]
        
    for {set i 1} {$i <= $number_state_var} {incr i} {
        set statevar$i [ttk::label $frstate.lb$i \
                -text [= "Initial state variable solid $i"]]
        set statevarval$i [ttk::entry $frstate.en$i \
                -textvariable "[$w give_uservar statevarval$i]"]      
        grid $frstate.lb$i -padx 2 -pady 2 -row $i -column 0
        grid $frstate.en$i -padx 2 -pady 2 -row $i -column 1   
    }
    
    $frstatevar.c create window 0 0 -anchor nw -window $frstate
    $frstatevar.c configure -scrollregion [$frstatevar.c bbox all]
    
    $f.framemat.note add $frstatevar -text [= "Initial state variables solid"] -underline 0 -padding 2
    
    grid $f -sticky news -ipadx 2 -ipady 2 -padx 2 -pady 2 
    grid columnconf $f 1 -weight 1
     
    grid columnconfigure $f 0 -weight 1
    grid rowconfigure $f 1 -weight 1
    
    set containerNode [$domNode selectNodes {container}] 
    
    for {set i 1} {$i <= $number_mat_params} {incr i} {
        if { $i<= 9 } {
            set ii "0$i" 
        } else {
            set ii "$i"
        }
        if { $containerNode == "" } {
            $w set_uservar_value matval$i 0.0         
        } else {             
            set n "material_parameter_solid_${ii}_"        
            set xp [format_xpath {value[@n=%s]} $n]
            set valueNode [$containerNode selectNodes $xp]
            if { $valueNode == "" } {
                $w set_uservar_value matval$i 0.0
            } else {
                $w set_uservar_value matval$i [$valueNode @v]            
            }           
        }
    }
    for {set i 1} {$i <= $number_state_var} {incr i} {
        if { $i<= 9 } {
            set ii "0$i" 
        } else {
            set ii "$i"
        }
        if { $containerNode == "" } {           
            $w set_uservar_value statevarval$i 0.0
        } else {                        
            set n "initial_state_variable_solid_${ii}_"
            set xp [format_xpath {value[@n=%s]} $n]
            set valueNode [$containerNode selectNodes $xp]
            if { $valueNode == "" } {
                $w set_uservar_value statevarval$i 0.0
            } else {
                $w set_uservar_value statevarval$i [$valueNode @v]  
            }
        }
    }
    

    bind $w <Return> +[list $w invokeok]
    bind $w <Escape> +[list MaterialParameters::invoke_escape $w]   
    
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
                MaterialParameters::save_data $w $domNode $number_mat_params $number_state_var
                catch { destroy $w }
                return
            }
        }        
    }
}

proc MaterialParameters::save_data { w domNode number_mat_params number_state_var } {   
    
    set doc $gid_groups_conds::doc      
         
    foreach node [list blockdata container value edit_command dependencies] {
        dom createNodeCmd element $node
    }
    set matList ""
    set statevarList ""
    for {set i 1} {$i <= $number_mat_params} {incr i} {
        set matval [$w give_uservar_value matval$i]                      
        set containerNode [$domNode selectNodes {container}]
        if { $containerNode == "" } {
            $domNode appendChildTag container \
                [list attributes() n "list_material_parameter_state_var" \
                    pn "List material parameters and state variables" state hidden]
            set containerNode [$domNode selectNodes {container}]
        } 
        if { $i<= 9 } {
            set ii "0$i" 
        } else {
            set ii "$i"
        }
        set isfound 0
        foreach iNode [$containerNode selectNodes {value}] {
            set n [$iNode @n]
            if { $n eq "material_parameter_solid_${ii}_"} {
                $iNode setAttribute v $matval
                lappend matList $n           
                set isfound 1     
            }                                    
        }
        if { !$isfound } {   
            set n "material_parameter_solid_${ii}_"
            set pn [lrange [split $n "_"] 0 end-1]        
            lappend matList $n
            $containerNode appendChildTag value \
                [list attributes() n $n pn $pn v $matval state hidden type number_mat_params]                    
        }                           
    }  
    for {set i 1} {$i <= $number_state_var} {incr i} {             
        set statevarval [$w give_uservar_value statevarval$i]
        
        set containerNode [$domNode selectNodes {container}]
        if { $containerNode == "" } {
            $domNode appendChildTag container \
                [list attributes() n "list_material_parameter_state_var" \
                    pn "List material parameters and state variables" state hidden]
            set containerNode [$domNode selectNodes {container}]
        } 
        if { $i<= 9 } {
            set ii "0$i" 
        } else {
            set ii "$i"
        }
        set isfound 0
        foreach iNode [$containerNode selectNodes {value}] {
            set n [$iNode @n]             
            if { $n eq "initial_state_variable_solid_${ii}_"} {
                $iNode setAttribute v $statevarval
                lappend statevarList $n  
                set isfound 1                
            }                         
        }
        if { !$isfound } {              
            set n "initial_state_variable_solid_${ii}_"
            set pn [lrange [split $n "_"] 0 end-1] 
            lappend statevarList $n        
            $containerNode appendChildTag value \
                [list attributes() n $n pn $pn v $statevarval state hidden type number_state_var]          
        }                           
    }
    foreach iNode [$containerNode selectNodes {value}] {
        set n [$iNode @n]
        set isfound 0
        if { $n in $matList || $n in $statevarList } {
            set isfound 1
        }      
        if {!$isfound} {
            $iNode delete
        }
    }
}