proc rmsd_slab_resting {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*2]
             set y1 [expr $highlim*5]

             set sel [atomselect top "name RM and x < $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open rmsdresting$no.csv w]
             
             for {set f $start} {$f < [molinfo top get numframes]} {incr f} {
                 set var $f
                 foreach index $index1 {
                     set ref2 [atomselect top "index $index" frame $start]
                     set sel2 [atomselect top "index $index" frame $f]
                     set rmsd [measure rmsd $sel2 $ref2]
                     append var ", " "$rmsd"
                     
                     unset ref2 
                     unset sel2 
                     unset rmsd 
                 }
             puts $outfile "$var"
             unset var     
             }
             close $outfile
         }
     }
}


proc slab_density_resting {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highlim]
            set x2 [expr $highlim*2]
            set x3 [expr $highlim*3]
            set x4 [expr $highlim*4]
            set x5 [expr $highlim*5]
            set x6 [expr $highlim*6]
            set x7 [expr $highlim*7]
            set y1 [expr $highlim*5]
            

            set outfile [open resting_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                
                set sel1 [atomselect top "name RM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2 [atomselect top "name RM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3 [atomselect top "name RM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4 [atomselect top "name RM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5 [atomselect top "name RM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6 [atomselect top "name RM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7 [atomselect top "name RM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                set num4 [$sel4 num]
                set num5 [$sel5 num]
                set num6 [$sel6 num]
                set num7 [$sel7 num]
                
                puts $outfile "$f, $num1, $num2, $num3, $num4, $num5, $num6, $num7"
                
                unset num1
                unset num2
                unset num3
                unset num4
                unset num5
                unset num6
                unset num7
                
                unset sel1
                unset sel2
                unset sel3
                unset sel4
                unset sel5
                unset sel6
                unset sel7
            }
            close $outfile
        }
        
    }
}

proc rmsd_slab_activated {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*2]
             set y1 [expr $highlim*5]

             set sel [atomselect top "name AM and x < $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open rmsdactivated$no.csv w]
             
             for {set f $start} {$f < [molinfo top get numframes]} {incr f} {
                 set var $f
                 foreach index $index1 {
                     set ref2 [atomselect top "index $index" frame $start]
                     set sel2 [atomselect top "index $index" frame $f]
                     set rmsd [measure rmsd $sel2 $ref2]
                     append var ", " "$rmsd"
                     
                     unset ref2 
                     unset sel2 
                     unset rmsd 
                 }
             puts $outfile "$var"
             unset var     
             }
             close $outfile
         }
     }
}



proc slab_density_activated {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highlim]
            set x2 [expr $highlim*2]
            set x3 [expr $highlim*3]
            set x4 [expr $highlim*4]
            set x5 [expr $highlim*5]
            set x6 [expr $highlim*6]
            set x7 [expr $highlim*7]
            set y1 [expr $highlim*5]
            

            set outfile [open activated_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                
                set sel1 [atomselect top "name AM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2 [atomselect top "name AM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3 [atomselect top "name AM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4 [atomselect top "name AM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5 [atomselect top "name AM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6 [atomselect top "name AM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7 [atomselect top "name AM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                set num4 [$sel4 num]
                set num5 [$sel5 num]
                set num6 [$sel6 num]
                set num7 [$sel7 num]
                
                puts $outfile "$f, $num1, $num2, $num3, $num4, $num5, $num6, $num7"
                
                unset num1
                unset num2
                unset num3
                unset num4
                unset num5
                unset num6
                unset num7
                
                unset sel1
                unset sel2
                unset sel3
                unset sel4
                unset sel5
                unset sel6
                unset sel7
            }
            close $outfile
        }
        
    }
}

proc slab_activated {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*5]
             set y1 [expr $highlim*5]

             set sel [atomselect top "name AM and x > $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open rmsdactivated$no.csv w]
             
             for {set f $start} {$f < [molinfo top get numframes]} {incr f} {
                 set var $f
                 foreach index $index1 {
                     set ref2 [atomselect top "index $index" frame $start]
                     set sel2 [atomselect top "index $index" frame $f]
                     set rmsd [measure rmsd $sel2 $ref2]
                     append var ", " "$rmsd"
                     
                     unset ref2 
                     unset sel2 
                     unset rmsd 
                 }
             puts $outfile "$var"
             unset var     
             }
             close $outfile
         }
     }
     
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highlim]
            set x2 [expr $highlim*2]
            set x3 [expr $highlim*3]
            set x4 [expr $highlim*4]
            set x5 [expr $highlim*5]
            set x6 [expr $highlim*6]
            set x7 [expr $highlim*7]
            set x8 [expr $highlim*8]
            set x9 [expr $highlim*9]
            set y1 [expr $highlim*5]
            

            set outfile [open activated_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                
                set sel1 [atomselect top "name AM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2 [atomselect top "name AM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3 [atomselect top "name AM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4 [atomselect top "name AM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5 [atomselect top "name AM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6 [atomselect top "name AM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7 [atomselect top "name AM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                set sel8 [atomselect top "name AM and x < $x8 and x > $x7 and y > 0 and y < $y1" frame $f]
                set sel9 [atomselect top "name AM and x < $x9 and x > $x8 and y > 0 and y < $y1" frame $f]
                set sel0 [atomselect top "name AM and x < $x9 and y > 0 and y < $y1" frame $f]
                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                set num4 [$sel4 num]
                set num5 [$sel5 num]
                set num6 [$sel6 num]
                set num7 [$sel7 num]
                set num8 [$sel8 num]
                set num9 [$sel9 num]
                set num0 [$sel0 num]
                
                puts $outfile "$f, $num1, $num2, $num3, $num4, $num5, $num6, $num7, $num8, $num9, $num0"
                
                unset num1
                unset num2
                unset num3
                unset num4
                unset num5
                unset num6
                unset num7
                unset num8
                unset num9
                unset num0
                
                unset sel1
                unset sel2
                unset sel3
                unset sel4
                unset sel5
                unset sel6
                unset sel7
                unset sel8
                unset sel9
                unset sel0
            }
            close $outfile
        }
        
    }
}

proc slab_resting {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*5]
             set y1 [expr $highlim*5]

             set sel [atomselect top "name RM and x < $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open rmsdresting$no.csv w]
             
             for {set f $start} {$f < [molinfo top get numframes]} {incr f} {
                 set var $f
                 foreach index $index1 {
                     set ref2 [atomselect top "index $index" frame $start]
                     set sel2 [atomselect top "index $index" frame $f]
                     set rmsd [measure rmsd $sel2 $ref2]
                     append var ", " "$rmsd"
                     
                     unset ref2 
                     unset sel2 
                     unset rmsd 
                 }
             puts $outfile "$var"
             unset var     
             }
             close $outfile
         }
     }
     
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highlim]
            set x2 [expr $highlim*2]
            set x3 [expr $highlim*3]
            set x4 [expr $highlim*4]
            set x5 [expr $highlim*5]
            set x6 [expr $highlim*6]
            set x7 [expr $highlim*7]
            set x8 [expr $highlim*8]
            set x9 [expr $highlim*9]
            set y1 [expr $highlim*5]
            

            set outfile [open resting_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                
                set sel1 [atomselect top "name RM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2 [atomselect top "name RM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3 [atomselect top "name RM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4 [atomselect top "name RM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5 [atomselect top "name RM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6 [atomselect top "name RM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7 [atomselect top "name RM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                set sel8 [atomselect top "name RM and x < $x8 and x > $x7 and y > 0 and y < $y1" frame $f]
                set sel9 [atomselect top "name RM and x < $x9 and x > $x8 and y > 0 and y < $y1" frame $f]
                set sel0 [atomselect top "name RM and x < $x9 and y > 0 and y < $y1" frame $f]
                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                set num4 [$sel4 num]
                set num5 [$sel5 num]
                set num6 [$sel6 num]
                set num7 [$sel7 num]
                set num8 [$sel8 num]
                set num9 [$sel9 num]
                set num0 [$sel0 num]
                
                puts $outfile "$f, $num1, $num2, $num3, $num4, $num5, $num6, $num7, $num8, $num9, $num0"
                
                unset num1
                unset num2
                unset num3
                unset num4
                unset num5
                unset num6
                unset num7
                unset num8
                unset num9
                unset num0
                
                unset sel1
                unset sel2
                unset sel3
                unset sel4
                unset sel5
                unset sel6
                unset sel7
                unset sel8
                unset sel9
                unset sel0
            }
            close $outfile
        }
        
    }
}

######################################################################################################
proc slab {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*2]
             set x3 [expr $highlim*3]
             set x4 [expr $highlim*4]
             set x5 [expr $highlim*5]
             set x6 [expr $highlim*6]
             set x7 [expr $highlim*7]
             set x8 [expr $highlim*8]
             set x9 [expr $highlim*9]
             set x10 [expr $highlim*10]
             set y1 [expr $highlim*5]

             set sela [atomselect top "name AM and x > $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1a [ lsort -integer -index 0 -increasing -unique [$sela get {index}]]
             set outfile1 [open rmsdactivated$no.csv w]
             set outfile2 [open activated_density$no.csv w]
             
             set selr [atomselect top "name RM and x < $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1r [ lsort -integer -index 0 -increasing -unique [$selr get {index}]]
             set outfile3 [open rmsdresting$no.csv w]
             set outfile4 [open resting_density$no.csv w]
             
             
             for {set f $start} {$f < [molinfo top get numframes]} {incr f} {
                 set vara $f
                 set varr $f
                 foreach index $index1a {
                    set refa [atomselect top "index $index" frame $start]
                    set selnewa [atomselect top "index $index" frame $f]
                    set rmsda [measure rmsd $selnewa $refa]
                    append vara ", " "$rmsda"
                     
                    unset refa 
                    unset selnewa
                    unset rmsda 
                    
                    }
                puts $outfile1 "$vara"
                unset vara
                
                #############################################################################################
                foreach index $index1r {
                    set refr [atomselect top "index $index" frame $start]
                    set selnewr [atomselect top "index $index" frame $f]
                    set rmsdr [measure rmsd $selnewr $refr]
                    append varr ", " "$rmsdr"
                     
                    unset refr 
                    unset selnewr 
                    unset rmsdr 
                    
                    }
                puts $outfile3 "$varr"
                unset varr  
                    
                set sel1a [atomselect top "name AM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2a [atomselect top "name AM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3a [atomselect top "name AM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4a [atomselect top "name AM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5a [atomselect top "name AM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6a [atomselect top "name AM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7a [atomselect top "name AM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                set sel8a [atomselect top "name AM and x < $x8 and x > $x7 and y > 0 and y < $y1" frame $f]
                set sel9a [atomselect top "name AM and x < $x9 and x > $x8 and y > 0 and y < $y1" frame $f]
                set sel10a [atomselect top "name AM and x < $x10 and x > $x9 and y > 0 and y < $y1" frame $f]
                set sel0a [atomselect top "name AM and x < $x10 and y > 0 and y < $y1" frame $f]
                     
                set num1a [$sel1a num]
                set num2a [$sel2a num]
                set num3a [$sel3a num]
                set num4a [$sel4a num]
                set num5a [$sel5a num]
                set num6a [$sel6a num]
                set num7a [$sel7a num]
                set num8a [$sel8a num]
                set num9a [$sel9a num]
                set num10a [$sel10a num]
                set num0a [$sel0a num] 
                    
                puts $outfile2 "$f, $num1a, $num2a, $num3a, $num4a, $num5a, $num6a, $num7a, $num8a, $num9a, $num10a, $num0a"
                
                unset num1a
                unset num2a
                unset num3a
                unset num4a
                unset num5a
                unset num6a
                unset num7a
                unset num8a
                unset num9a
                unset num10a
                unset num0a
                
                unset sel1a
                unset sel2a
                unset sel3a
                unset sel4a
                unset sel5a
                unset sel6a
                unset sel7a
                unset sel8a
                unset sel9a
                unset sel10a
                unset sel0a
                ############################################################################################################
                set sel1r [atomselect top "name RM and x < $x1 and x > 0   and y > 0 and y < $y1" frame $f]
                set sel2r [atomselect top "name RM and x < $x2 and x > $x1 and y > 0 and y < $y1" frame $f]
                set sel3r [atomselect top "name RM and x < $x3 and x > $x2 and y > 0 and y < $y1" frame $f]
                set sel4r [atomselect top "name RM and x < $x4 and x > $x3 and y > 0 and y < $y1" frame $f]
                set sel5r [atomselect top "name RM and x < $x5 and x > $x4 and y > 0 and y < $y1" frame $f]
                set sel6r [atomselect top "name RM and x < $x6 and x > $x5 and y > 0 and y < $y1" frame $f]
                set sel7r [atomselect top "name RM and x < $x7 and x > $x6 and y > 0 and y < $y1" frame $f]
                set sel8r [atomselect top "name RM and x < $x8 and x > $x7 and y > 0 and y < $y1" frame $f]
                set sel9r [atomselect top "name RM and x < $x9 and x > $x8 and y > 0 and y < $y1" frame $f]
                set sel10r [atomselect top "name RM and x < $x10 and x > $x9 and y > 0 and y < $y1" frame $f]
                set sel0r [atomselect top "name RM and x < $x10 and y > 0 and y < $y1" frame $f]
                     
                set num1r [$sel1r num]
                set num2r [$sel2r num]
                set num3r [$sel3r num]
                set num4r [$sel4r num]
                set num5r [$sel5r num]
                set num6r [$sel6r num]
                set num7r [$sel7r num]
                set num8r [$sel8r num]
                set num9r [$sel9r num]
                set num10r [$sel10r num]
                set num0r [$sel0r num] 
                    
                puts $outfile4 "$f, $num1r, $num2r, $num3r, $num4r, $num5r, $num6r, $num7r, $num8r, $num9r, $num10r, $num0r"
                
                unset num1r
                unset num2r
                unset num3r
                unset num4r
                unset num5r
                unset num6r
                unset num7r
                unset num8r
                unset num9r
                unset num10r
                unset num0r
                
                unset sel1r
                unset sel2r
                unset sel3r
                unset sel4r
                unset sel5r
                unset sel6r
                unset sel7r
                unset sel8r
                unset sel9r
                unset sel10r
                unset sel0r
                  
            }
            close $outfile1
            close $outfile2
            close $outfile3
            close $outfile4
        }
        
    }
}
