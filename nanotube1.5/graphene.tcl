# Generate graphene sheets
#
# $Id: graphene.tcl,v 1.9 2013/04/15 15:27:10 johns Exp $
#

proc ::Nanotube::graphene_usage { } {
    vmdcon -info "Usage: graphene -lx <length> -ly <length> -type <edge type> \[-nlayers <number of layers>\]  \[-b <0|1>\] \[-a <0|1>\] \[-d <0|1>\] \[-i <0|1>\] \[-noh <0.0-1.0>\] \[-cc <blength>\] \[-ma <C-C|B-N>\]"
    vmdcon -info "  <length> is the edge length in nanometers"
    vmdcon -info "  <edge type> is the type of edge (armchair or zigzag)"
    vmdcon -info "  <blength> is the length of a bond in nanometers (default: 0.1418)"    
    vmdcon -info "  <number of layers> is the number of layers of graphene (default: 1)"
    vmdcon -info "  -ma C-C/B-N selects between carbon and boron-nitride (default: C-C)"
    vmdcon -info "  -b 0/1 turns generation of bonds off/on (default: on)"
    vmdcon -info "  -a 0/1 turns generation of angles off/on (default: on)"
    vmdcon -info "  -d 0/1 turns generation of dihedrals off/on (default: on)"
    vmdcon -info "  -i 0/1 turns generation of impropers off/on (default: on)"
    vmdcon -info "  -nbo 0.0-100.0 add the percentage of -O- groups to GRA (default: 0)"
    vmdcon -info "  -noh 0.0-100.0 add the percentage of -OH groups to GRA (default: 0)"
    vmdcon -info "  -ncooh 0.0-100.0 add the percentage of -COOH groups to GRA (default: 0)"
    vmdcon -info "  -pbcBO 0/1 turns pbc condition off/on for -O- (default: no)"
    vmdcon -info "  -pbcOH 0/1 turns pbc condition off/on for -OH (default: no)"
    vmdcon -info "  -pbcCOOH 0/1 turns pbc condition off/on for -COOH(default: no)"
    vmdcon -info "  -addH  0/1  adding all edge hydrogens(default: no)"
    vmdcon -info "  -savemol2  0/1  save molecules to gar.mol2(default: no)"
    vmdcon -info "  The -a/-d/-i flags only have an effect if -b 1 is used"
}

# a simple shuffle algorithm: https://wiki.tcl-lang.org/page/Shuffle+a+list
proc ::Nanotube::shuffle {list} {
    set n [llength $list]
    for {set i 1} {$i < $n} {incr i} {
        set j [expr {int(rand() * $n)}]
        set temp [lindex $list $i]
        lset list $i [lindex $list $j]
        lset list $j $temp
    }
    return $list
}

# rotate residue in rad
proc ::Nanotube::rotateFunc {sel angle} {
    set cood1 [lindex [$sel get {x y z}] 0]
    set coord2 [lindex [$sel get {x y z}] 1]
    set offset [vecsub {0 0 0} $cood1]
    # move C-group to zero
    $sel move [transoffset $offset]
    # random rotate C-group along V12
    set V12 [vecsub $cood1 $coord2]
    $sel move [trans axis $V12 $angle rad]
    # move back C-group to orgin
    set offset [vecsub $cood1 {0 0 0}]
    $sel move [transoffset $offset]
}

# rotate residue along fix axis
proc ::Nanotube::rotateAxis {sel angle axis} {
    set cood1 [lindex [$sel get {x y z}] 0]
    set offset [vecsub {0 0 0} $cood1]
    # move C-group to zero
    $sel move [transoffset $offset]
    # rotate along $axis
    $sel move [transabout $axis $angle rad]
    # move back C-group to orgin
    set offset [vecsub $cood1 {0 0 0}]
    $sel move [transoffset $offset]
}

proc ::Nanotube::graphene_core { args } {
    # Check if proper #arguments was given
    set n_args [llength $args]
    if { [expr fmod($n_args,2)] } {
        vmdcon -err "graphene: wrong number of arguments"
        vmdcon -err ""
        graphene_usage 
        return -1
    }
    if { ($n_args < 8) || ($n_args > 40) } {
        vmdcon -err "graphene: wrong number of arguments"
        vmdcon -err ""
        graphene_usage
        return -1
    }

    # build a full topology by default
    set cmdline(-b) 1
    set cmdline(-a) 1
    set cmdline(-d) 1
    set cmdline(-i) 1 
    set cmdline(-nbo) 0         ; #the percentage of -O- 
    set cmdline(-pbcBO) 0       ; #if use pbc for -O-
    set cmdline(-noh) 0         ; #the percentage of -OH 
    set cmdline(-pbcOH) 0       ; #if use pbc for -OH
    set cmdline(-ncooh) 0       ; #the percentage of -COOH 
    set cmdline(-pbcCOOH) 0     ; #if use pbc for -COOH
    set cmdline(-addH) 0        ; # add edge Hydrogens
    set cmdline(-savemol2) 0    ; # save mol2
    set cmdline(-ma) C-C
    set cmdline(-cc) 0.1418
    set cmdline(-space) 3.35

    for { set i 0} {$i < $n_args} {incr i 2} {
        set key [lindex $args $i]
        set val [lindex $args [expr $i + 1]]
        set cmdline($key) $val
    }

    # Check if mandatory options are defined
    foreach a {-lx -ly -type} {
        if { ![info exists cmdline($a)] } {
            vmdcon -err "graphene: required flag '$a' is missing\n"
            graphene_usage
            return -1
        }
    }
  
    if { [info exists cmdline(-nlayers)] } {
        set nlayers $cmdline(-nlayers)
    } else {
        set nlayers 1
    }
    # Set graphene parameters
    set lx $cmdline(-lx)
    set ly $cmdline(-ly)
    set type $cmdline(-type)
    set a [expr {10.0*$cmdline(-cc)}]
    set pi 3.14159265358979323846
    set space $cmdline(-space)
    set nbo $cmdline(-nbo)
    set p_bo $nbo
    set pbcBO $cmdline(-pbcBO)
    set noh $cmdline(-noh)
    set p_oh $noh
    set pbcOH $cmdline(-pbcOH)
    set ncooh $cmdline(-ncooh)
    set p_cooh $ncooh
    set pbcCOOH $cmdline(-pbcCOOH)
    set addH $cmdline(-addH)
    set savemol2 $cmdline(-savemol2)

    #Check that input is reasonable
    if { $lx <=0 || $ly <= 0 } {
        vmdcon -err "graphene: Edge length must be positive"
        return -1
    }
    if { ($type != "armchair") && ($type != "zigzag")} {
        vmdcon -err "graphene: Type must be either 'armchair' or 'zigzag'"
        return -1
    }

    if {($cmdline(-ma) != {C-C}) && ($cmdline(-ma) != {B-N})} {
        vmdcon -err "graphene: Material (-ma) must be either 'C-C' or 'B-N'"
        return -1
    }
    
    # check input parameters for functional groups
    if {$noh+$ncooh>100.0} {
        error "graphene: the sum percentage of all groups must be <= 100.0"
        return -1
    }

    #Number of unit cells
    if {$type=="armchair"} {
        set Lx_cell [expr 2*$a*sin(60*$pi/180)]
        set Ly_cell [expr 3*$a]
    } else {
        set Lx_cell [expr 3*$a]
        set Ly_cell [expr 2*$a*sin(60*$pi/180)]
    }

    set Nx_cell [expr ceil($lx*10/$Lx_cell)]
    set Ny_cell [expr ceil($ly*10/$Ly_cell)]

    #Index min/max
    set i 0

    #Generate unit cell coordinates
    if {$type=="armchair"} {
        set r1 "0 0 0"
        set r2 "[expr -$a*sin(60*$pi/180)] [expr $a*cos(60*$pi/180)] 0"
        set r3 [vecadd $r2 "0 $a 0"]
        set r4 "0 [expr 2*$a] 0"
        set l_shift "0 $a 0"
    } else {
        set r1 "0 0 0"
        set r2 "[expr -$a*cos(60*$pi/180)] [expr $a*sin(60*$pi/180)] 0"
        set r3 "$a 0 0"
        set r4 [vecadd $r2 "[expr 2*$a] 0 0"]
        set l_shift "[expr $a/2] [expr $a*sin(60*$pi/180)] 0" 
    }

    #Generate graphene coordinates
    set xyzlist {}
    set num_atoms 0
    #edge C idx list
    set edgelist {}
    for {set k 0} { $k < $nlayers} {incr k} {
        for {set j 0} { $j < $Ny_cell } {incr j} {
            for {set i 0} { $i < $Nx_cell } {incr i} {
                set r_shift "[expr $i*$Lx_cell] [expr $j*$Ly_cell] [expr $k*$space]"
                if {[expr $k%2]!=0} {set r_shift [vecadd $r_shift $l_shift]}
                lappend xyzlist [vecadd $r1 $r_shift] 
                lappend xyzlist [vecadd $r2 $r_shift]
                lappend xyzlist [vecadd $r3 $r_shift]
                lappend xyzlist [vecadd $r4 $r_shift]
                # get edge atoms index list
                if {$j==0} {
                    lappend edgelist [expr $num_atoms]
                    if {$type!="armchair"} { lappend edgelist [expr $num_atoms+2] }
                } 
                if {$j==$Ny_cell-1} {
                    lappend edgelist [expr $num_atoms+3]
                    if {$type!="armchair"} { lappend edgelist [expr $num_atoms+1]}
                } 
                if {$i==0} {
                    lappend edgelist [expr $num_atoms+1]
                    if {$type=="armchair"} { lappend edgelist [expr $num_atoms+2]} 
                }
                if {$i==$Nx_cell-1} {
                    lappend edgelist [expr $num_atoms+3] 
                    if {$type=="armchair"} { lappend edgelist [expr $num_atoms] } 
                }
                incr num_atoms 4
            }
        }
    }
    set edgelist [lsort -unique $edgelist]
    set old_atoms $num_atoms
    set name    [lrepeat $old_atoms C]
    set resname [lrepeat $old_atoms GRA]
    set segid   [lrepeat $old_atoms SHT]
    set element [lrepeat $old_atoms C]
#    set type    [lrepeat $old_atoms CA]
    set type    [lrepeat $old_atoms C]
    set mass    [lrepeat $old_atoms 12.0107]
    set radius  [lrepeat $old_atoms 1.7]
    set chain   [lrepeat $old_atoms X]
    set charge  [lrepeat $old_atoms 0.0]
    
    set blist {}
    if {$nbo>0} {
        # first build pure graphene to get the bonds only for adding -O- group
        vmdcon -info "-------------------------------------------------------------"
        set mol [mol new atoms $old_atoms]
        animate dup $mol
        set sel [atomselect $mol all]
        $sel set {x y z} $xyzlist
        foreach key {name resname segid element type mass radius chain charge} \
                value "[list $name $resname $segid $element $type $mass $radius $chain $charge]" {
            $sel set $key $value
        }
        mol bondsrecalc $mol
        mol reanalyze $mol
        set g_bonds [$sel getbonds]
        set g_nbonds [$sel get numbonds]
        $sel delete
        mol delete $mol
        
        # get the number of bonds for the graphene
        for {set i 0} {$i<[llength $g_nbonds]} {incr i} {
            if {[lindex $g_nbonds $i]>0} {
                foreach idx [lindex $g_bonds $i] {
                    lappend blist [lsort -integer [list $i $idx]]
                }
            }
        }
        set blist [lsort -unique $blist]
        set nbonds  [llength $blist]
        # check it, if the number of -O- > the number of bonds
        if {$nbo/100.0 * $old_atoms > $nbonds} {
            error "Too many -O-, at least <= $nbonds"
        }
        vmdcon -info "Total bonds= $nbonds"
        vmdcon -info "-------------------------------------------------------------"
    }
    
    if {$noh+$ncooh+$nbo>0} {
        # the number of -O- to int
        set nbo [expr int($nbo/100.0 * $old_atoms)]
        # the number of -OH to int
        set noh [expr int($noh/100.0 * $old_atoms)]
        # the number of -COOH to int
        set ncooh [expr int($ncooh/100.0 * $old_atoms)]
        vmdcon -info "Need add no= $nbo, noh= $noh, ncooh= $ncooh"
        # generate -OH atoms coordinates
        set idxlist [list]
        # roughly check enough C to add -OH, -COOH, -O-(one -O- occupy two C atoms)
        set nsum [expr $noh+$ncooh+$nbo*2]
        set nfree $old_atoms
        set NoPBC [expr {!$pbcOH && !$pbcCOOH && !$pbcBO}]
        if {$NoPBC} {incr nfree -[llength $edgelist]}
        if {$nsum > $nfree} {
            set nfixed 100.0
            if {$NoPBC} {
                set nfixed [expr 100.0 * ($old_atoms-[llength $edgelist]) / $old_atoms]
            }
            error "The sum of percentage(%) must be <= [format "%.1f" $nfixed]"
        }
        # 0-based graphene C atom index for add groups
        for {set i 0} {$i < $old_atoms} {incr i} {
            # no pbc condition
            if {$NoPBC} {
                if {$i ni $edgelist} {lappend idxlist $i}
            } else {
                lappend idxlist $i
            }
        }

        set idxlist [shuffle $idxlist]
        set idxC_O [list]  ;# {C1-O-C2 order}
        # add -O-
        if {$nbo>0} {
            set newblist $blist
            # exclude edege bonds part from newblist
            if {!$pbcBO} {
                set newblist {}
                foreach bond $blist {
                    set ia [lindex $bond 0]
                    set ib [lindex $bond 1]
                    if {$ia in $idxlist && $ib in $idxlist} {
                        lappend newblist [list $ia $ib]
                    }
                }
            }
            # unorder all bonds
            set newblist [shuffle $newblist]
            # put -O- to bonds upper or lower
            set count 0
            set nb [llength $newblist]
            for {set i 0} {$i < $nb} {incr i} {
                if {$count==$nbo} { break }
                
                set bond [lindex $newblist $i]
                set ia [lindex $bond 0]
                set ib [lindex $bond 1]
                # center point of C-C
                set cent [ vecscale 0.5 [vecadd [lindex $xyzlist $ia] [lindex $xyzlist $ib]] ]
                set coordO [vecadd $cent {0 0 1.2}]
                if {[expr rand()]>0.5} {
                    set coordO [vecadd $cent {0 0 -1.2}]
                }
                
                # remove has added two C atoms of -O- and return new list
                set idxa [lsearch $idxlist $ia]
                set idxb [lsearch $idxlist $ib]
                if {$idxa >= 0 && $idxb >= 0} {
                    set idxlist [lreplace $idxlist $idxa $idxa]
                    # must re-search it because idxlist has changed
                    set idxb [lsearch $idxlist $ib]
                    set idxlist [lreplace $idxlist $idxb $idxb]
                    lappend xyzlist $coordO
                    lappend name OC
                    lappend idxC_O [list $ia [expr $old_atoms+$count] $ib]
                    incr count
                }
                
                # check again...
                if {$i==$nb-1 && $count < $nbo} {
                    error "Not found enough C-C to add -O-"
                }
            }
            
            incr num_atoms $nbo
        }
        
        # add -OH
        set idxC_OH [list]
        set oldlist $idxlist
        set count 0
        for {set i 0} {$i < [llength $oldlist]} {incr i} {
            # has got noh enough
            if {$count==$noh} { break }
            # skip the index
            if {!$pbcOH && [lindex $oldlist $i] in $edgelist} {
                continue
            }
        
            set coordOH {{0 0 0} {0 0 1.47} {0 0.91 1.77}}
            set natoms [llength $coordOH]
            # -coordOH
            if {[expr rand()]>0.5} {
                for {set j 0} {$j < $natoms} {incr j} {
                    lset coordOH $j [vecscale -1.0 [lindex $coordOH $j]]
                }
            }
            set targC [lindex $xyzlist [lindex $oldlist $i]]
            set vec [vecsub $targC [lindex $coordOH 0]]
            for {set j 1} {$j < $natoms} {incr j} {
                lappend xyzlist [vecadd $vec [lindex $coordOH $j]]
            }
            lappend name OA HA
            
            # incr count 
            incr count
            # save C index
            lappend idxC_OH [lindex $oldlist $i]
            # remove has added -OH  C atom and return new list
            set idx [lsearch $idxlist [lindex $oldlist $i]]
            if {$idx < 0} {
                error "remove is error"
            }
            set idxlist [lreplace $idxlist $idx $idx]
            
            # check again...
            if {$i==[llength $oldlist]-1 && $count < $noh} {
                error "Not found enough C to add -OH"
            }
        }
        # update atoms
        incr num_atoms [expr $noh * 2]
        
        # add COOH
        set idxC_COOH [list]
        set oldlist $idxlist
        set count 0
        for {set i 0} {$i < [llength $oldlist]} {incr i} {
            # has got ncooh enough
            if {$count==$ncooh} { break }
            # skip the index
            if {!$pbcCOOH && [lindex $oldlist $i] in $edgelist} {
                continue
            }
        
            set coordCOOH {{0 0 0} {0 0 1.56} {0 0.99 2.23} {0 -1.23 2.06} {0 -1.21 3.03}}
            set natoms [llength $coordCOOH]
            # -coordCOOH
            if {[expr rand()]>0.5} {
                for {set j 0} {$j < $natoms} {incr j} {
                    lset coordCOOH $j [vecscale -1.0 [lindex $coordCOOH $j]]
                }
            }
            set targC [lindex $xyzlist [lindex $oldlist $i]]
            set vec [vecsub $targC [lindex $coordCOOH 0]]
            for {set j 1} {$j < $natoms} {incr j} {
                lappend xyzlist [vecadd $vec [lindex $coordCOOH $j]]
            }
            # COOH
            lappend name CB OB OB HB
            
            # incr count 
            incr count
            # save C indx 
            lappend idxC_COOH [lindex $oldlist $i]
            # remove has added -COOH C atom and return new list
            set idx [lsearch $idxlist [lindex $oldlist $i]]
            if {$idx < 0} {
                error "remove is error"
            }
            set idxlist [lreplace $idxlist $idx $idx]
            
            # check again...
            if {$i==[llength $oldlist]-1 && $count < $ncooh} {
                error "Not found enough C to add -COOH"
            }
        }
        # update atoms
        incr num_atoms [expr $ncooh * 4]
        
        # add ...
    }
    
    # add Edege Hydrogens
    if {$addH} {
        # fake H atoms
        set nadd [expr [llength $edgelist]+2]
        for {set i 0} {$i < $nadd} {incr i} {
            lappend xyzlist {1000 1000 1000}
            lappend name HH
        }
        incr num_atoms $nadd
    }

    #Create new molecule with one frame
    set mol [mol new atoms $num_atoms]
    animate dup $mol
    set sel [atomselect $mol all]
    set asel [atomselect $mol {index % 2 == 0}]
    set bsel [atomselect $mol {index % 2 == 1}]
    set mat Graphene
    # set box
    set box "[expr $Nx_cell*$Lx_cell] [expr $Ny_cell*$Ly_cell] [expr $nlayers*$space]"
    vmdcon -info "CRYST1= $box"
    molinfo $mol set {a b c} $box
    
    #Set default values for all atoms
    if {$cmdline(-ma) == {C-C}} {
        # append functional group info
        for {set i $old_atoms} {$i < $num_atoms} {incr i} {
            set atom_name [string index [lindex $name $i] 0]
            set atom_res GRA
            set atom_segid SHT
            set atom_element C
#            set atom_type CA
            set atom_type C
            set atom_mass 12.0107
            set atom_radius 1.7
            set atom_chain X
            set atom_charge 0.10
            if {$atom_name=="O"} {
                set atom_element O
                set atom_type O
                set atom_mass 16.0
#                set atom_radius 1.52
                set atom_radius 0.0  ; # set up zero turns off auto-rebonds
            } elseif { $atom_name=="H"} {
                set atom_element H
                set atom_type H
                set atom_mass 1.008
#                set atom_radius 1.0
                set atom_radius 0.0 ; # set up zero turns off auto-rebonds
            }
            lappend resname $atom_res
            lappend segid   $atom_segid
            lappend element $atom_element
            lappend type    $atom_type
            lappend mass    $atom_mass
            lappend radius  $atom_radius
            lappend chain   $atom_chain
            lappend charge  $atom_charge
        }
        foreach key {name resname segid element type mass radius chain charge} \
                value "[list $name $resname $segid $element $type $mass $radius $chain $charge]" {
            $sel set $key $value
        }
    } elseif {$cmdline(-ma) == {B-N}} {
        set mat {Boron Nitride}
        foreach key {name resname segid element type mass radius chain charge} value {B BNS SHT B B 10.811 1.7265 X 1.05} {
            $asel set $key $value
        }
        foreach key {name resname segid element type mass radius chain charge} value {N BNS SHT N N 14.0067 1.6825 X -1.05} {
            $bsel set $key $value
        }
    }
    $asel delete
    $bsel delete

    $sel set {x y z} $xyzlist


    #Add representation for molecule
    if {$type=="armchair"} {mol rename $mol "Armchair $mat Sheet"}
    if {$type=="zigzag"} {mol rename $mol "Zigzag $mat Sheet"}

    # only build topology information that is enabled
    if {($cmdline(-b) == "on") || ($cmdline(-b) == 1)} {
        # stash away current radius information for reliable bond searching
        set oldrad [$sel get radius]
        # is need??
#        $sel set radius [expr {1.2*$a}]
        mol bondsrecalc $mol
        
        # add conect bonds to selected pure C atoms
        set bonds [$sel getbonds]
        if {$nbo+$noh+$ncooh>0} {
            # add -O- two bonds
            foreach atoms $idxC_O {
                set idxC1 [lindex $atoms 0]
                set idxO  [lindex $atoms 1]
                set idxC2 [lindex $atoms 2]
                # C1-O
                set newbond [lindex $bonds $idxC1]
                lappend newbond $idxO
                lset bonds $idxC1 $newbond
                # C2-O
                set newbond [lindex $bonds $idxC2]
                lappend newbond $idxO
                lset bonds $idxC2 $newbond
                # O-C1, O-C2
                lset bonds $idxO [list $idxC1 $idxC2]
            }
            # update old atoms because add -OH
            incr old_atoms $nbo
            
            # add -OH bonds
            for {set i 0} {$i<$noh} {incr i} {
                # C-O bonds
                set idxC [lindex $idxC_OH $i]
                set idxO [expr $old_atoms+2*$i]
                set idxH [expr $old_atoms+2*$i+1]
                set newbond [lindex $bonds $idxC]
                lappend newbond $idxO
                lset bonds $idxC $newbond
                # O-C, O-H
                set newbond [list $idxC $idxH]
                lset bonds $idxO $newbond
                # H-O bonds
                lset bonds $idxH $idxO 
                
                # rotate -OH
                set OH [atomselect $mol "index $idxC $idxO $idxH"]
                rotateFunc $OH [expr rand()*2*$pi]
            }
            # update old atoms because add -OH
            incr old_atoms [expr $noh * 2]
            
            # add -COOH bonds
            for {set i 0} {$i<$ncooh} {incr i} {
                set idxC  [lindex $idxC_COOH $i] ;# origin GRA C atom
                set idxCA [expr $old_atoms+4*$i] ; # C atom of -COOH
                set idxO1 [expr $old_atoms+4*$i+1]
                set idxO2 [expr $old_atoms+4*$i+2]
                set idxH  [expr $old_atoms+4*$i+3]
                # C-CA, CA-C, CA-O1, CA-O2
                set newbond [lindex $bonds $idxC] ; # get old bonds
                lappend newbond $idxCA
                lset bonds $idxC $newbond
                lset bonds $idxCA [list $idxC $idxO1 $idxO2]
                # O1-CA, O2-CA, O2-H    
                lset bonds $idxO1 $idxCA  
                lset bonds $idxO2 [list $idxCA $idxH]
                # H-O2
                lset bonds $idxH $idxO2
                
                # rotate -COOH
                set COOH [atomselect $mol "index $idxC $idxCA $idxO1 $idxO2 $idxH"]
                rotateFunc $COOH [expr rand()*2*$pi]
            }
            # update old atoms because add -COOH
            incr old_atoms [expr $ncooh * 4]
        }
        
        # move edge hydrogens to target positions and add bonds
        if {$addH} {
            # get updated coordinates
            set xyzlist [$sel get {x y z}]
            # edge C atom index
            foreach idx $edgelist {
                set bond [lindex $bonds $idx]
                set numbond [llength $bond]
                set atom0xyz [lindex $xyzlist $idx]
                # the C need add One H
                if {$numbond==2} {
                    set atom1xyz [lindex $xyzlist [lindex $bond 0]]
                    set atom2xyz [lindex $xyzlist [lindex $bond 1]]
                    set V10 [vecsub $atom0xyz $atom1xyz]
                    set V20 [vecsub $atom0xyz $atom2xyz]
                    set Vnorm [vecadd $V10 $V20]
                    set Vnorm [vecscale [expr 1.08 * 1.0/[veclength $Vnorm]] $Vnorm]
                    # change Hydroge to positon
                    lset xyzlist $old_atoms [vecadd $Vnorm $atom0xyz]

                    # add C-H bonds, H-C bonds
                    set newbond [lindex $bonds $idx]
                    lappend newbond $old_atoms
                    lset bonds $idx $newbond
                    lset bonds $old_atoms $idx
                    
                    # update atoms
                    incr old_atoms
                } else {   ;# the C need add two H
                    # TODO: more faster, relpace $sel get and set operation
                    set atom1xyz [lindex $xyzlist [lindex $bond 0]]
                    set V10 [vecsub $atom0xyz $atom1xyz]
                    set Vnorm [vecscale [expr 1.08 * 1.0/[veclength $V10]] $V10]
                    for {set i -1} {$i < 2} {incr i 2} {
                        lset xyzlist $old_atoms [vecadd $Vnorm $atom0xyz]
                        $sel set {x y z} $xyzlist
                        # rotate pos to get two C-H bonds along Z axis
                        set CH [atomselect $mol "index $idx $old_atoms"]
                        rotateAxis $CH [expr $i*$pi/3.0] {0 0 1}
                        # must update xyzlist ...
                        set xyzlist [$sel get {x y z}]
                        
                        # add C-H, H-C bonds
                        set newbond [lindex $bonds $idx]
                        lappend newbond $old_atoms
                        lset bonds $idx $newbond
                        lset bonds $old_atoms $idx
                        
                        incr old_atoms
                    }
                }
            }
            $sel set {x y z} $xyzlist
        }
        
        # added new bonds
        $sel setbonds $bonds
        
        # restore original radius
        $sel set radius $oldrad
        
        # set bond types. this will also trigger the flag that
        # the bonds will be written out through molfile.
        ::TopoTools::retypebonds $sel

        if {($cmdline(-a) == "on") || ($cmdline(-a) == 1)} {
            ::TopoTools::guessangles $sel
        }
        if {($cmdline(-d) == "on") || ($cmdline(-a) == 1)} {
            ::TopoTools::guessdihedrals $sel
        }
        if {($cmdline(-i) == "on") || ($cmdline(-a) == 1)} {
            ::TopoTools::guessimpropers $sel {tolerance 5}
        }
    }
    mol reanalyze $mol
    ::TopoTools::adddefaultrep $mol

    $sel set resid [$sel get serial]
    # save mol2?
    if {$savemol2} { 
        set systemTime [clock seconds]
        set t [clock format $systemTime -format %Y%m%d_%H-%M-%S]
        $sel writemol2 "GRA_$p_oh%OH_$p_cooh%COOH_${t}.mol2" 
    }
    
    $sel delete
    return $mol
}

# insert the textmode command variant into the default namespace
interp alias {} graphene {} ::Nanotube::graphene_core
