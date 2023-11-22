# Carbon Nanostructure Builder GUI
#
# $Id: gui.tcl,v 1.7 2013/04/15 15:27:10 johns Exp $
#

package require topotools 1.3

namespace eval ::Nanotube:: {
    variable w            ; # handle of the GUI's toplevel widget
    variable version  1.5 ; # version number of this plugin

    variable ccbond 0.1418; # length of bonds in nanometers
    variable ccorbn C-C   ; # choose between Carbon and BN system

    variable l   5        ; # nanotube length in nanometers
    variable n   5        ; # nanotube chirality parameter n
    variable m  10        ; # nanotube chirality parameter m

    variable lx  5        ; # x-length of graphene sheet in nanometers
    variable ly  5        ; # y-length of graphene sheet in nanometers
    variable lz  3.0      ; # distance between layers in stack
    variable type armchair; # graphene sheet boundary type.
    variable nlayers    1 ; # number of sheets in stack
    variable space     3.35 ; # space of sheets in stack, Angstrom
    variable stacking  AB ; # AA, AB, ABC

    variable bonds      1 ; # generate bonds  
    variable angles     1 ; # generate angles
    variable dihedrals  1 ; # generate dihedrals
    variable impropers  1 ; # generate impropers
    
    # for functional group
    variable noh  0  ; # the precent of -OH
    variable pbcOH  0 ; # if use pbc condition for add -OH (default: no/0)
    variable ncooh  0
    variable pbcCOOH 0
    variable nbo      0
    variable pbcBO    0
    variable addH 0
    variable savemol2 0
}
package provide nanotube $::Nanotube::version

proc ::Nanotube::nanotube_gui {} {
    variable w
    variable version

    variable ccbond
    variable ccorbn

    variable l
    variable n
    variable m

    variable lx
    variable ly
    variable lz
    variable type
    variable nlayers
    variable space
    variable stacking

    variable bonds
    variable angles
    variable dihedrals
    variable impropers
    
    # for functional group
    variable noh
    variable pbcOH
    variable ncooh
    variable pbcCOOH
    variable nbo
    variable pbcBO
    variable addH
    variable savemol2

    if { [winfo exists .nanotube] } {
        wm deiconify $w
        return
    }

    set w [toplevel ".nanotube"]
    wm title $w "Carbon Nanostructure Builder"
    wm resizable $w no no
    set row 0

    #Add a menubar
    frame $w.menubar -relief raised -bd 2 -padx 10
    grid  $w.menubar -padx 1 -column 0 -columnspan 4 -row $row -sticky ew
    menubutton $w.menubar.help -text "Help" -underline 0 \
        -menu $w.menubar.help.menu
    $w.menubar.help config -width 5
    pack $w.menubar.help -side right
    menu $w.menubar.help.menu -tearoff no
    $w.menubar.help.menu add command -label "About" \
        -command {tk_messageBox -type ok -title "About Carbon Nanostructure Builder" \
                      -message "Tool for building selected\ncarbon nanostructures.\n\nVersion $::Nanotube::version\n\n(c) 2009-2011 \nby Robert R. Johnson\n <robertjo@physics.upenn.edu>\nand\nAxel Kohlmeyer\n<akohlmey@gmail.com>"}
    $w.menubar.help.menu add command -label "Help..." \
        -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/nanotube"
    incr row

    grid [label $w.topoopts -justify center -relief raised -text "Topology Building Options:"] \
        -row $row -column 0 -columnspan 4 -sticky nsew
    incr row
    grid [label $w.matlabel -text "Material: "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [menubutton $w.material -width 7 -relief raised -bd 2 \
              -direction flush -text "test text" -menu $w.material.menu \
              -textvariable ::Nanotube::ccorbn] \
        -row $row -column 3 -columnspan 1 -sticky ew
    menu $w.material.menu -tearoff no
    $w.material.menu add radiobutton -value C-C -command {set ::Nanotube::ccbond 0.1418} -label C-C -variable ::Nanotube::ccorbn
    $w.material.menu add radiobutton -value B-N -command {set ::Nanotube::ccbond 0.1446} -label B-N -variable ::Nanotube::ccorbn
    incr row
    grid [label $w.cclabel -text "Length of bond (nm): "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.ccbond -width 7 -textvariable ::Nanotube::ccbond] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row
    grid [checkbutton $w.topobonds -variable ::Nanotube::bonds -text Bonds \
              -command [namespace code \
                            {
                                if {$bonds} {
                                    $w.topoangles configure -state active
                                    $w.topodihedrals configure -state active
                                    $w.topoimpropers configure -state active
                                } {
                                    $w.topoangles configure -state disabled
                                    $w.topodihedrals configure -state disabled
                                    $w.topoimpropers configure -state disabled}
                            } ] ] \
        -row $row -column 0 -sticky nsew
    grid [checkbutton $w.topoangles -variable ::Nanotube::angles -text Angles] \
        -row $row -column 1 -sticky nsew
    grid [checkbutton $w.topodihedrals -variable ::Nanotube::dihedrals -text Dihedrals] \
        -row $row -column 2 -sticky nsew
    grid [checkbutton $w.topoimpropers -variable ::Nanotube::impropers -text Impropers] \
        -row $row -column 3 -sticky nsew
 
    incr row
    grid [label $w.tubeopts -justify center -relief raised -text "Nanotube Building Options:"] \
        -row $row -column 0 -columnspan 4 -sticky nsew
    incr row
    grid [label $w.nlabel -text "Nanotube chiral index n: "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.n -width 7 -textvariable ::Nanotube::n] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [label $w.mlabel -text "Nanotube chiral index m: "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.m -width 7 -textvariable ::Nanotube::m] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [label $w.llabel -text "Nanotube length (nm): "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.l -width 7 -textvariable ::Nanotube::l] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [button $w.gotube -text "Generate Nanotube" -foreground "blue" \
              -command [namespace code \
                            { vmdcon -info \
                                  "calling nanotube_core -l $l -n $n -m $m -b $bonds -a $angles -d $dihedrals -i $impropers -cc $ccbond -ma $ccorbn"
                              nanotube_core -l "$l" -n "$n" -m "$m" -b "$bonds" \
                                  -a "$angles" -d "$dihedrals" -i "$impropers" -cc "$ccbond" -ma "$ccorbn"
                            } ]] -row $row -column 0 -columnspan 4 -sticky nsew
    incr row

    grid [label $w.sheetopts -justify center -relief raised -text "Graphene Sheet Building Options:"] \
        -row $row -column 0 -columnspan 4 -sticky nsew
    incr row
    grid [label $w.lxlabel -text "Edge length along x (nm): "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.lx -width 7 -textvariable ::Nanotube::lx] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [label $w.lylabel -text "Edge length along y (nm): "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.ly -width 7 -textvariable ::Nanotube::ly] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [label $w.layerlabel -text "Number of layers: "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.nlayers -width 7 -textvariable ::Nanotube::nlayers] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row
    
    grid [label $w.layerspacelabel -text "Distance of layers (Angstrom): "] \
        -row $row -column 0 -columnspan 3 -sticky w
    grid [entry $w.space -width 7 -textvariable ::Nanotube::space] \
        -row $row -column 3 -columnspan 1 -sticky ew
    incr row

    grid [label $w.typelabel -text "Graphene edge type: "] \
        -row $row -column 0 -columnspan 2 -sticky w
    grid [radiobutton $w.armchair -text "Armchair" -value "armchair" -variable ::Nanotube::type] \
        -row $row -column 2 -columnspan 1 -sticky w
    grid [radiobutton $w.zigzag -text "Zigzag"   -value "zigzag"   -variable ::Nanotube::type] \
        -row $row -column 3 -columnspan 1 -sticky w
        
    # functional group menu
    incr row 
    grid [label $w.nolabel -text "-O- (%): "] -row $row -column 0 -columnspan 1 -sticky w
    grid [entry $w.nbo -width 7 -textvariable ::Nanotube::nbo] \
        -row $row -column 1 -columnspan 1 -sticky w
    grid [checkbutton $w.pbcBO -variable ::Nanotube::pbcBO -text "use PBC"] \
        -row $row -column 2 -sticky w
        
    incr row 
    grid [label $w.nohlabel -text "-OH (%): "] -row $row -column 0 -columnspan 1 -sticky w
    grid [entry $w.noh -width 7 -textvariable ::Nanotube::noh] \
        -row $row -column 1 -columnspan 1 -sticky w
    grid [checkbutton $w.pbcOH -variable ::Nanotube::pbcOH -text "use PBC"] \
        -row $row -column 2 -sticky w
        
    incr row 
    grid [label $w.ncoohlabel -text "-COOH (%): "] -row $row -column 0 -columnspan 1 -sticky w
    grid [entry $w.ncooh -width 7 -textvariable ::Nanotube::ncooh] \
        -row $row -column 1 -columnspan 1 -sticky w
    grid [checkbutton $w.pbcCOOH -variable ::Nanotube::pbcCOOH -text "use PBC"] \
        -row $row -column 2 -sticky w
    incr row 
    
    grid [checkbutton $w.savemol2 -variable ::Nanotube::savemol2 -text "save to GRAxxx.mol2"] \
        -row $row -column 0 -sticky w
    grid [checkbutton $w.addH -variable ::Nanotube::addH -text "add Edge Hydrogen" \
        -command [namespace code \
                 {
                    if {$addH} {
                        set pbcBO 0
                        set pbcOH 0
                        set pbcCOOH 0
                        $w.pbcBO configure -state disabled
                        $w.pbcOH configure -state disabled
                        $w.pbcCOOH configure -state disabled
                    } else {
                        $w.pbcBO configure -state active
                        $w.pbcOH configure -state active
                        $w.pbcCOOH configure -state active
                    }
                 }
                 ] ]\
        -row $row -column 2 -sticky nsew
        
    incr row
    grid [button $w.gosheet -text "Generate Sheet(s) " -foreground "blue" \
              -command [namespace code \
                            { vmdcon -info \
                                  "calling graphene_core -lx $lx -ly $ly -type $type \
                                  -nlayers $nlayers -b $bonds -a $angles -d $dihedrals \
                                  -i $impropers -noh $noh -pbcOH $pbcOH \
                                  -ncooh $ncooh -pbcCOOH $pbcCOOH -nbo $nbo -pbcBO $pbcBO \
                                  -addH $addH -savemol2 $savemol2 -cc $ccbond -space $space \
                                  -ma $ccorbn"
                              graphene_core -lx "$lx" -ly "$ly" -type "$type" \
                                  -nlayers "$nlayers" -b "$bonds" -a "$angles" \
                                  -d "$dihedrals" -i "$impropers" -noh "$noh" \
                                  -pbcOH $pbcOH -ncooh $ncooh -pbcCOOH $pbcCOOH \
                                  -nbo $nbo -pbcBO $pbcBO -space $space \
                                  -addH $addH -savemol2 $savemol2 -cc "$ccbond" -ma "$ccorbn"
                            } ]] -row $row -column 0 -columnspan 4 -sticky nsew
}

proc nanotube_tk {} {
  ::Nanotube::nanotube_gui
  return $::Nanotube::w
}
