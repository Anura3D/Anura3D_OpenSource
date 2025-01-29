#######################################################
# STAGES MANAGER WINDOW 
#######################################################

proc Anura3D::StagesManager { {widget ""} } {
    global Priv
    
    package require fulltktree

    set doc $gid_groups_conds::doc
    set root [$doc documentElement]           

    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    }      
    
    Anura3D::CreateImages
           
    set cols {
        { 14 { Status } left text 0 }
        { 15 { Stage ID }  left text 0 }  
        { 100 { Calculation state } left text 0 }
    }   
    
    package require img::png
    if {$widget != ""} {
        set w $widget
    } else {
        set w .gid.t
    }
    catch {destroy $w}
    toplevel $w      
    wm title $w [= "Stages manager"]
    
    set t [fulltktree $w.a -columns $cols -width 585 -height 350 -showlines 0 -indent 0 -selectmode extended \
            -sensitive_cols all -showbuttons 0 -buttonpress_open_close 0 -returnhandler 0]
             
    set xp {container[@n='stages']/blockdata[@n='stage']}  
       
    foreach stageselectedNode [$root selectNodes $xp]  { 
        set calcstate [Anura3D::ReadGeneralInfo "calculation_state" $stageselectedNode]       
        set prevcalcstate [Anura3D::ReadGeneralInfo "prevcalcstate" $stageselectedNode]               
        if {$calcstate eq "wrongcalc" || $calcstate eq "shouldnotcalc"} {
#             $calcstate eq "OKcalc"
            set calcstatenew "notcalc"                     
            Anura3D::UpdateGeneralInfo "calculation_state" $calcstatenew $stageselectedNode 
            Anura3D::UpdateGeneralInfo "prevcalcstate" $calcstate $stageselectedNode 
            Anura3D::UpdateGeneralInfo "mark_calculate" 1 $stageselectedNode                 
        }
    } 
    
    set ipos 0     
    set message [= "Click on the status icon to change its current status."]
    foreach stagesselectedNode [$root selectNodes $xp]  {  
        set calcstate [Anura3D::ReadGeneralInfo "calculation_state" $stagesselectedNode]             
        set stagename [$stagesselectedNode @name]
        set prevcalcstate [Anura3D::ReadGeneralInfo "prevcalcstate" $stagesselectedNode]
        
        if {$calcstate eq "notcalc"} {
            set image "navforward16"
            set help [= "Stage should be calculated."]
            if {$prevcalcstate == "notcalc"} {
                set imageold "forward-16"
                set helpold [= "Stage should NOT be calculated."]
            } elseif {$prevcalcstate eq "OKcalc"} {
                set imageold "actcheck16"
                set helpold [= "Stage has been calculated successfully."]
            } elseif {$prevcalcstate eq "wrongcalc"} {
                set imageold "actcross16"
                set helpold [= "Calculation of this stage failed."]
            } elseif {$prevcalcstate eq "shouldnotcalc"} {
                set imageold "forward-16"
                set helpold [= "Stage should NOT be calculated."]
            }
        } elseif {$calcstate eq "OKcalc"} {
            set image "actcheck16"
            set help [= "Stage has been calculated successfully."]
            if {$prevcalcstate == "notcalc"} {
                set imageold "navforward16"
                set helpold [= "Stage should be calculated."]
            } elseif {$prevcalcstate eq "OKcalc"} {
                set imageold "navforward16"
                set helpold [= "Stage should be calculated"]
            } elseif {$prevcalcstate eq "wrongcalc"} {
                set imageold "navforward16"
                set helpold [= "Stage should be calculated."]              
            } elseif {$prevcalcstate eq "shouldnotcalc"} {
                set imageold "navforward16"
                set helpold [= "Stage should be calculated."]            
            }           
        } elseif {$calcstate eq "wrongcalc"} {
            set image "actcross16"
            set help [= "Calculation of this stage failed."]
            set imageold "navforward16"
            set helpold [= "Stage should be calculated."]
        } elseif {$calcstate eq "shouldnotcalc"} {
            set image "forward-16"
            set help [= "Stage should NOT be calculated."]
            set imageold "navforward16"
            set helpold [= "Stage should be calculated."]
        }
                      
        set item [$w.a insert end [list "" $stagename $help $helpold $imageold]] 
        
        menubutton $w.a.m$ipos -image $image -menu $w.a.m$ipos.m \
            -relief flat -bd 1 -highlightthickness 0 \
            -background white -foreground blue -activeforeground red -width 16
        menu $w.a.m$ipos.m -tearoff 0
        set DialogWin::user($ipos,status) 1
        set DialogWin::user($ipos,image) $image
        
        $w.a.m$ipos.m add command -label "Change status" -command \
            [list change $w.a $item $image $imageold $help $helpold $ipos $stagename]          
        
        gid_groups_conds::register_popup_help $w.a.m$ipos $help            
        
        $w.a item style set $item 0 window
        $w.a item style set $item 1 text
        $w.a item style set $item 2 text
        $w.a item element configure $item 0 e_window -window $w.a.m$ipos
        $w.a configure -itemheight 22
        incr ipos          
    }  
        
    ttk::label $w.note -wraplength 500 -text \
        [= "Note: Click on the status icons to change the calculation state for each stage."]   
    
    ttk::button $w.calc -text [= "Calculate"] -command [list Anura3D::ok $w $ipos 1]
    ttk::button $w.c -text [= "Cancel"] -command [list Anura3D::cancel $w]
    set modelname [GiD_Info Project Modelname]    
    
    grid $t - - - -sticky nsew -padx 3 -pady 3
    grid $w.note $w.calc $w.c -sticky w -padx 3 -pady 3
    grid conf $w.calc $w.c -sticky ew         
    
    grid rowconfigure $w 0 -weight 1
    grid rowconfigure $w 1 -weight 0   
    grab $w
    focus $w
    update
}

proc change { tree itemList image imageold help helpold idx stagename } { 
    variable doc
    
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]  
    
    set stageList ""
    set xp {container[@n='stages']/blockdata[@n='stage']}      
  
    foreach item $itemList {
        set mb [$tree item element cget $item 0 e_window -window]
        if {$DialogWin::user($idx,status)} { 
            $mb configure -image $imageold
            gid_groups_conds::register_popup_help $mb $helpold
            set DialogWin::user($idx,status) 0 
            set DialogWin::user($idx,image) $imageold 
            if {$image ne "actcheck16"} { continue }
            set jdx 0; set jitem 1        
        } else {
            $mb configure -image $image
            gid_groups_conds::register_popup_help $mb $help
            set DialogWin::user($idx,status) 1  
            set DialogWin::user($idx,image) $image
        }
    }    
}

proc Anura3D::ok { w ipos callcalculate } {     
    variable doc
          
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]  
    
    set ipos 0
    set stagecalcList ""
    set xp {container[@n='stages']/blockdata[@n='stage']}     
    foreach stageselectedNode [$root selectNodes $xp]  {    
        set stagename [$stageselectedNode @name]   
        set markcalculate [Anura3D::ReadGeneralInfo "mark_calculate" $stageselectedNode]                     
        set calcstate [Anura3D::ReadGeneralInfo "calculation_state" $stageselectedNode]   
        set prevcalcstate [Anura3D::ReadGeneralInfo "prevcalcstate" $stageselectedNode]                                               
                
        if {$DialogWin::user($ipos,image) eq "navforward16"} {
            set DialogWin::user($ipos,markcalc) 1    
        } else {
            set DialogWin::user($ipos,markcalc) 0    
        }                                                                      
        if {$DialogWin::user($ipos,markcalc) == "1"} {                                                                             
            set markcalc [Anura3D::ReadGeneralInfo "mark_calculate" $stageselectedNode]  
            set calcstate [Anura3D::ReadGeneralInfo "calculation_state" $stageselectedNode]                                                                                    
            lappend stagecalcList $stagename
        }
        incr ipos      
    }                                     
    set ipos 0
    foreach stageNode [$root selectNodes $xp]  {         
        set stagename [$stageNode @name]   
        Anura3D::UpdateGeneralInfo "mark_calculate" $DialogWin::user($ipos,markcalc) $stageNode            
        set calcstate [Anura3D::ReadGeneralInfo "calculation_state" $stageNode]   
        if {$DialogWin::user($ipos,status)} { 
            set calcstatenew $calcstate
        } else {
            if {$calcstate eq "notcalc"} {
                set calcstatenew "shouldnotcalc"
            } elseif {$calcstate eq "OKcalc" || $calcstate eq "wrongcalc" || $calcstate eq "shouldnotcalc"} {
                set calcstatenew "notcalc"
            } 
        }
        Anura3D::UpdateGeneralInfo "calculation_state" $calcstatenew $stageNode 
        Anura3D::UpdateGeneralInfo "prevcalcstate" $calcstate $stageNode                               
        incr ipos            
    }    
    
    gid_groups_conds::actualize_conditions_window    
    
    grab .gid
    focus .gid  
    
    if { $callcalculate } {
        if { $stagecalcList == ""} {
            snit_messageBox -parent .gid -message [= "A stage should be marked to calculate."] 
            return 
        } 
        destroy $w                                      
        Anura3D::Calculate $stagecalcList 
    } else {
        destroy $w
    }      
}

proc Anura3D::cancel { w } {
    destroy $w
    grab .gid
    focus .gid
    update  
}

proc Anura3D::CreateImages {} {
    image create photo navforward16 -data {
        R0lGODlhEAAQAIUAAPwCBAwyTBRObAw2VDR+nCRKZOzy/KTe7Pz+/KTK3Nzu
        /Lze7FS+1AyexAyuzBSavAyOtBSmzOTy/BRqjNTm9IzO5ETS3ETa5By61Ayi
        xByixBRmjAQGDBxCXGSivCySrCSWtBTC3AQOHAQWHAxWdEze7AQKFBRCXAwq
        PAQCBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAAAZj
        QIBwSCwahYGjUjBQGgWEpHNYMBCaT4G2UDggos+EwmBYMBpf6VBgYDgeEMgj
        IpmoAQVKxXLBPDIXGhscRB0eHyAgDSGBGyJFASMiIiMkJYImUwAnmJqbjp4A
        KCmhAKSlTn5BACH+aENyZWF0ZWQgYnkgQk1QVG9HSUYgUHJvIHZlcnNpb24g
        Mi41DQqpIERldmVsQ29yIDE5OTcsMTk5OC4gQWxsIHJpZ2h0cyByZXNlcnZl
        ZC4NCmh0dHA6Ly93d3cuZGV2ZWxjb3IuY29tADs=
    }
    
    image create photo actcross16 -data {
        R0lGODlhEAAQAIIAAASC/PwCBMQCBEQCBIQCBAAAAAAAAAAAACH5BAEAAAAA
        LAAAAAAQABAAAAMuCLrc/hCGFyYLQjQsquLDQ2ScEEJjZkYfyQKlJa2j7AQn
        MM7NfucLze1FLD78CQAh/mhDcmVhdGVkIGJ5IEJNUFRvR0lGIFBybyB2ZXJz
        aW9uIDIuNQ0KqSBEZXZlbENvciAxOTk3LDE5OTguIEFsbCByaWdodHMgcmVz
        ZXJ2ZWQuDQpodHRwOi8vd3d3LmRldmVsY29yLmNvbQA7
    }
    
    image create photo actcheck16 -data {
        R0lGODlhEAAQAIIAAPwCBMT+xATCBASCBARCBAQCBEQCBAAAACH5BAEAAAAA
        LAAAAAAQABAAAAM2CLrc/itAF8RkdVyVye4FpzUgJwijORCGUhDDOZbLG6Nd
        2xjwibIQ2y80sRGIl4IBuWk6Af4EACH+aENyZWF0ZWQgYnkgQk1QVG9HSUYg
        UHJvIHZlcnNpb24gMi41DQqpIERldmVsQ29yIDE5OTcsMTk5OC4gQWxsIHJp
        Z2h0cyByZXNlcnZlZC4NCmh0dHA6Ly93d3cuZGV2ZWxjb3IuY29tADs=
    }
    image create photo forward-16 -data {
        R0lGODlhEAAQAIUAAPwCBHy2hIS6jHy2jPz+/JzGlPT67HyyhISyhNzy3Oz2
        7Ozy5OTy3OTy5HSidJTOhFyWXKTSnLzetMzqzNTq1FyKXIzKhJzSlMTixNzq
        1Nzu3Hy6dER+VITCfITGfJzOlFSmZERmTDxyRERuTGSmXGyybHS2bIS+fJTG
        jDRiRDxeRHS+fER2XDxmVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAAAZq
        QIBwSCwah4GjUjBQGgOEpFFQAASYAgOhWTwQAIiEYsFoRIuOBwASkUwoDcWW
        WFFXLJELBpPRaKQAFRsAHB0eFh4FHyAhRCIQACMkJSYnKCApRSqZKSEhECuM
        TkMcmKNDLJmnQi2rrq9+QQAh/mhDcmVhdGVkIGJ5IEJNUFRvR0lGIFBybyB2
        ZXJzaW9uIDIuNQ0KqSBEZXZlbENvciAxOTk3LDE5OTguIEFsbCByaWdodHMg
        cmVzZXJ2ZWQuDQpodHRwOi8vd3d3LmRldmVsY29yLmNvbQA7
    } 

    image create photo back-22 -data {
        R0lGODlhFgAWAIUAAPwCBAw2VCRGZAxCZGyavExmjHyatOTy9CxihISevPz+
        /KzO3BRylAw+XAQCBDRWbPz6/FzC3CSuzDyexJzO5Mzq9CxSdAQOFISmxNzu
        9HTS5BSmxAyexDSuzJTa7Mzu9Kzi7GS21CRmjAQOHHSWtLze7AyWvHzG3BRi
        hAQKFCTO3BS+1AyixBSWvBSOtBSStAQWJBSixDzW5BTC3BSqzBS21CTC1ETW
        3AQSHEze7BRqlBRmjAQCDBR+pBRefBRSdCH5BAEAAAAALAAAAAAWABYAAAal
        QIBwSCwaj8ikMqBcMpvHgGAANQYIhWdVGDAcENQtIJBQLBgNx0MQaDuQXcgh
        IplQDhBIxXKJYiAZGhscHR4VHyAhIiNWJBklGhIbJoQnFCcTKIxFKSgbKiss
        Ji0mJi4vLiYoMEcXKDEyMzQ1Nje2NisoOEg4KDU5K6g6OwwoKAN9SCOeMmgw
        z884PEq9PT4NYkPLP9jZQikN3d4AKVrjKePp3gZBACH+aENyZWF0ZWQgYnkg
        Qk1QVG9HSUYgUHJvIHZlcnNpb24gMi41DQqpIERldmVsQ29yIDE5OTcsMTk5
        OC4gQWxsIHJpZ2h0cyByZXNlcnZlZC4NCmh0dHA6Ly93d3cuZGV2ZWxjb3Iu
        Y29tADs=
    }

    image create photo forward-22 -data {
        R0lGODlhFgAWAIUAAPwCBAw2VAQCBBxCXDR+nIS21Aw+XJTC1Nzu/KzO3Pz+
        /Nzq9Pz6/MTe7KTW5FzC1Nzu9CRKZMzi7IzK3Lzi7LTe7HzG3Gy+3AyuzAye
        xFzC3DRSbHy+1Dy61CSqzAySvAyStLze7IzO5AyGrETa5ByixBRmjCTC1ETS
        3BTC3Bx2nAyWvEze7AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAWABYAAAaY
        QIBwSCwaj8hkMqBsBgTN5IAAjRoDBaq1aDggtMuAWDzoJhTgY+CwYLgZDccD
        wkgXI5IJZVGxXDAZGnR2QxsLHB0PHRgeHyAZDyFfVUQDCyIgIyCPIB+QJCUm
        lEMBEiInKCQnKSkeKSQeomoqJrUmKiArKSwZsmoCwMEBGCyxo1EGHr3HUQEE
        vltCBtDRAAbMW0zV29xDBkEAIf5oQ3JlYXRlZCBieSBCTVBUb0dJRiBQcm8g
        dmVyc2lvbiAyLjUNCqkgRGV2ZWxDb3IgMTk5NywxOTk4LiBBbGwgcmlnaHRz
        IHJlc2VydmVkLg0KaHR0cDovL3d3dy5kZXZlbGNvci5jb20AOw==
    }  
  
    image create photo pencil-22 -data {
        R0lGODlhFgAWAIMAAASC/IQCBPwCBPyChMQCBPzCxAQCBPz+/MzKzISChKyq
        rDQyNEQCBAAAAAAAAAAAACH5BAEAAAAALAAAAAAWABYAAARYEMhJ6wxiEMtp
        IAWxddwXiqRlikSQeiAbuC+wirNR322gv7zcLobzDU+9XypoBBKTR1lz+RTW
        Dgip8nUwZK1XLyIx5XoVicX2RUAo1DVKi7GOBxjxfNwQAQAh/mhDcmVhdGVk
        IGJ5IEJNUFRvR0lGIFBybyB2ZXJzaW9uIDIuNQ0KqSBEZXZlbENvciAxOTk3
        LDE5OTguIEFsbCByaWdodHMgcmVzZXJ2ZWQuDQpodHRwOi8vd3d3LmRldmVs
        Y29yLmNvbQA7
    }  
    image create photo eraser-22 -data {
        R0lGODlhFgAWAIMAAPwCBMTC/AQChAQCBISC/PzCxPz+/MQCBIQCBAAAAAAA
        AAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAWABYAAARfEMhJKwhi2F0D
        IRnHeR+oiRRZmqikrmD7wgI6r+EI4+d2l7lNwWDYsTiFw4FI61mSymXxE3xG
        o8xqBXpVGrQUbveAcG7H0bKowEar12wx2SyMp+lI+7s1ie/5fX8tBhEAIf5o
        Q3JlYXRlZCBieSBCTVBUb0dJRiBQcm8gdmVyc2lvbiAyLjUNCqkgRGV2ZWxD
        b3IgMTk5NywxOTk4LiBBbGwgcmlnaHRzIHJlc2VydmVkLg0KaHR0cDovL3d3
        dy5kZXZlbGNvci5jb20AOw==
    }  
    image create photo up-22 -data {
        R0lGODlhFgAWAIUAAPwCBAw2VAQCBHSWtBRmjAQOHISmxNzu9BSmxBRihHya
        tPz6/Lze7CTO3BSixHTS5BTC3DzW5ByyzPz+/OTy9AyexEze7ByixGyavKzO
        3FzC3AyWvBS+1BR+pAQKFCRGZExmjCxihBRylCSuzBSWvBS21BSOtBRSdAw+
        XAxCZDyexDSyzCTC1JzO5JTa7DSuzETW3BRqlAQWJDRWbOT2/Mzq9HzG3JzS
        5Kzi7BSStGS21CxSdCRmjAQOFAQSHAAAACH5BAEAAAAALAAAAAAWABYAAAae
        QIBwSCwaj8ikcqkMCJjHwIBQgBIDhgMiUbUGFAtGw0GFfheHByQi4S6/E8pD
        UoFYLm5kAEPJaBAVGxIcER0JHlEfICEiIxUkGyUmIgknKIhXASkonCorgSwm
        KQGcKE9IAi0uLxUwMTJWMzQ1NiYwBLBQHws1N7avXgs4NjkcCblMATU6KhvG
        yG87PAnUKV1MAj0+2zIFp1bg4eJJdkEAIf5oQ3JlYXRlZCBieSBCTVBUb0dJ
        RiBQcm8gdmVyc2lvbiAyLjUNCqkgRGV2ZWxDb3IgMTk5NywxOTk4LiBBbGwg
        cmlnaHRzIHJlc2VydmVkLg0KaHR0cDovL3d3dy5kZXZlbGNvci5jb20AOw==
    }
    image create photo down-22 -data {
        R0lGODlhFgAWAIUAAPwCBAw2VCRKZDRSbBxCXJTC1Mzi7Nzq9NTm9Bx2nAQC
        BNzu9JzG3Hy+1HzG3IzO5BRmjPz6/LTe7Dy61AyStCTC1FzC1AyGrETS3ETC
        1ETa5BRulAyuzBRylAw+XMTe7Gy+3CSqzAyexBTC3DR+nIS21KTW5Nzu/KzO
        3FzC3Pz+/ByixEze7AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAWABYAAAaR
        QIBwSCwaj8ikcnkMBAQDgjPAFAYKhsMBkVBUAYEFo+F4QLzVQEQyoVTOX/XB
        csHA0+vMRbNBMwkRDhxuHX5GTlIeHh8gISIjFAEeiVRECiQlDAUmgxQjIhwi
        JHdFlycoKSIUFCEjGiGkRpcqCxYijxorsUezcxYsuoZJsxLAu0qXB7DCTJfH
        VQrMX9PU1Uh0QQAh/mhDcmVhdGVkIGJ5IEJNUFRvR0lGIFBybyB2ZXJzaW9u
        IDIuNQ0KqSBEZXZlbENvciAxOTk3LDE5OTguIEFsbCByaWdodHMgcmVzZXJ2
        ZWQuDQpodHRwOi8vd3d3LmRldmVsY29yLmNvbQA7
    }

    image create photo accelerate -data {
        R0lGODlhFgAWAIUAAPwCBDQyNFxeXAQCBMTGxOzm7CwqLLy2vPTy9Pz+/Ly6
        vCQiJLSytLS2tLSutOTi5MzGzKSepIyKjJSOlKSmpMzKzJyanIyOjBwaHIyG
        jISGhJSSlISChBQSFJyenIR+hGxubDw+PHRydHR2dEQ+RHx6fERCRAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAWABYAAAa6
        QIBwSCwaj8hkIIBcJgEBweAYnTYJUmMAa9USClniFtwlGg6IRFhoUKTXwwWj
        0FB3F46Hwl6UQyISfAB+EROCQgsUFRYSF3yJEIyBaxgWDBkaGRtclQwSHBIb
        EGEdGx4fGhcOICEDGBsWHBmqIq1CHRIWGRMMIyRTHRy6Er22tyONq8YdJRe0
        xkIDwr2/QwMfliMmZQADIxasZd4e4UYDIr7c59rc0eVFA+/m0EQD9PDt0flP
        /P3+BkEAACH+aENyZWF0ZWQgYnkgQk1QVG9HSUYgUHJvIHZlcnNpb24gMi41
        DQqpIERldmVsQ29yIDE5OTcsMTk5OC4gQWxsIHJpZ2h0cyByZXNlcnZlZC4N
        Cmh0dHA6Ly93d3cuZGV2ZWxjb3IuY29tADs=
    }
}

