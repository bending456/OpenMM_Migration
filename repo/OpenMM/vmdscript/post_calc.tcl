proc postcalc_slab { num_sim_series repeat minFrame maxFrame highLim} {
    # prefix : prefix number
    # num_sim_series: a number of simulation in a series
    # repeat: a number of repetition 
    for {set k 1 } {$k <= $num_sim_series} {incr k} {
        for {set i 0 } {$i < $repeat} {incr i} {
            set no [expr $k + $num_sim_series*$i]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highLim/3]
            set x2 [expr $highLim*2/3]

            set selr1 [atomselect top "name AR and x < $x1"]
            set indexr1 [lsort -integer -index 0 -increasing -unique [$selr1 get {index}]]

            set selr2 [atomselect top "name AR and x > $x1 and x < $x2"]
            set indexr2 [lsort -integer -index 0 -increasing -unique [$selr2 get {index}]]

            set selr3 [atomselect top "name AR and x > $x2"]
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

            set countr2 0
            set xr2 0
            set yr2 0 
            set rr2 0
            set outfile2 [open test_r2_$no.csv w]

            foreach index $indexr2 {

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

                set countr2 [expr $countr2 + 1]

                puts $outfile2 "$countr2, $x_dist, $y_dist, $r_dist"

                $sel delete

            }
            $selr2 delete
            close $outfile2

            set countr3 0
            set xr3 0
            set yr3 0 
            set rr3 0
            set outfile3 [open test_r3_$no.csv w]

            foreach index $indexr3 {

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

                set countr3 [expr $countr3 + 1]

                puts $outfile3 "$countr3, $x_dist, $y_dist, $r_dist"

                $sel delete

            }
            $selr3 delete
            close $outfile3

            set selr1 [atomselect top "name HR and x < $x1"]
            set indexr1 [lsort -integer -index 0 -increasing -unique [$selr1 get {index}]]

            set selr2 [atomselect top "name HR and x > $x1 and x < $x2"]
            set indexr2 [lsort -integer -index 0 -increasing -unique [$selr2 get {index}]]

            set selr3 [atomselect top "name HR and x > $x2"]
            set indexr3 [lsort -integer -index 0 -increasing -unique [$selr3 get {index}]]

            set countr1 0
            set xr1 0
            set yr1 0 
            set rr1 0
            set outfile1 [open test_a1_$no.csv w]

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

            set countr2 0
            set xr2 0
            set yr2 0 
            set rr2 0
            set outfile2 [open test_a2_$no.csv w]

            foreach index $indexr2 {

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

                set countr2 [expr $countr2 + 1]

                puts $outfile2 "$countr2, $x_dist, $y_dist, $r_dist"

                $sel delete

            }
            $selr2 delete
            close $outfile2

            set countr3 0
            set xr3 0
            set yr3 0 
            set rr3 0
            set outfile3 [open test_a3_$no.csv w]

            foreach index $indexr3 {

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

                set countr3 [expr $countr3 + 1]

                puts $outfile3 "$countr3, $x_dist, $y_dist, $r_dist"

                $sel delete

            }
            $selr3 delete
            close $outfile3   
        } 
    }         
}