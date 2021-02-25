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
