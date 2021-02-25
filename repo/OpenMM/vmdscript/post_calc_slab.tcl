proc postcalc_slab { num_sim_series repeat minFrame maxFrame highLim} {
    # prefix : prefix number
    # num_sim_series: a number of simulation in a series
    # repeat: a number of repetition 
    for {set k 1 } {$k <= $num_sim_series} {incr k} {
        for {set i 0 } {$i < $repeat} {incr i} {
            set no [expr $k + $num_sim_series*$i]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highLim]
            set x2 [expr $highLim*2]

            set selr1 [atomselect top "name RM and x < $x1" frame 1]
            set indexr1 [lsort -integer -index 0 -increasing -unique [$selr1 get {index}]]

            set selr2 [atomselect top "name RM and x > $x1 and x < $x2" frame 1]
            set indexr2 [lsort -integer -index 0 -increasing -unique [$selr2 get {index}]]

            set selr3 [atomselect top "name RM and x > $x2" frame 1]
            set indexr3 [lsort -integer -index 0 -increasing -unique [$selr3 get {index}]]

            set countr1 0
            set outfile1 [open test_r1_$no.csv w]

            foreach index $indexr1 {

                set sel [atomselect top "index $index"]
                $sel frame $minFrame
                set coordmin [measure center $sel]
                set xcoordmin [lindex $coordmin 0]
                set ycoordmin [lindex $coordmin 1]

                $sel frame $maxFrame
                set coordmax [measure center $sel]
                set xcoordmax [lindex $coordmax 0]
                set ycoordmax [lindex $coordmax 1]

                set x_dist [expr abs($xcoordmax - $xcoordmin)]
                set y_dist [expr abs($ycoordmax - $ycoordmin)]
                set r_dist [expr sqrt(($x_dist*$x_dist)+($y_dist*$y_dist))]

                set countr1 [expr $countr1 + 1]

                puts $outfile1 "$countr1, $x_dist, $y_dist, $r_dist"

                $sel delete

            }
            $selr1 delete
            close $outfile1

        } 
    }         
}

 proc postcalc {start highlim num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set x1 [expr $highlim]
             set x2 [expr $highlim*2]
             set y1 [expr $highlim*5]

             set sel [atomselect top "name RM and x < $x1 and x > 0 and y > 0 and y < $y1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open test$no.csv w]
             
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