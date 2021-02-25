proc postcalc { num_sim_series repeat minFrame maxFrame highLim} {
    # prefix : prefix number
    # num_sim_series: a number of simulation in a series
    # repeat: a number of repetition 
    for {set k 1 } {$k <= $num_sim_series} {incr k} {
        for {set i 0 } {$i < $repeat} {incr i} {
            set no [expr $k + $num_sim_series*$i]
            mol load pdb test$no.pdb dcd test$no.dcd
            set x1 [expr $highLim/2]

            set selr1 [atomselect top "name RM and x < $x1" frame 1]
            set indexr1 [lsort -integer -index 0 -increasing -unique [$selr1 get {index}]]

            set selr2 [atomselect top "name RM and x > $x1" frame 1]
            set indexr2 [lsort -integer -index 0 -increasing -unique [$selr2 get {index}]]

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

            set selr1 [atomselect top "name AM and x < $x1" frame 1]
            set indexr1 [lsort -integer -index 0 -increasing -unique [$selr1 get {index}]]

            set selr2 [atomselect top "name AM and x > $x1" frame 1]
            set indexr2 [lsort -integer -index 0 -increasing -unique [$selr2 get {index}]]

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

        }
    }
}

proc rmsd {num_sim_series repeat minFrame maxFrame highLim} {

     for {set k 1 } {$k <= $num_sim_series} {incr k} {
        for {set i 0 } {$i < $repeat} {incr i} {
        
            set no [expr $k + $num_sim_series*$i]
            mol load pdb test$no.pdb dcd test$no.dcd
            
            set x1 [expr $highLim/2]
                       
            set reference1 [atomselect top "name RM and x < $x1" frame 1]
            set comparison1 [atomselect top "name RM and x < $x1" frame 1]
            
            set reference2 [atomselect top "name RM and x > $x1" frame 1]
            set comparison2 [atomselect top "name RM and x > $x1" frame 1]
            
            #set reference3 [atomselect top "name AM and x < $x1" frame 1]
            #set comparison3 [atomselect top "name AM and x < $x1" frame 1]
            
            #set reference4 [atomselect top "name AM and x > $x1" frame 1]
            #set comparison4 [atomselect top "name AM and x > $x1" frame 1]

            set outfile1 [open rmsd_r1_$no.csv w]
            set outfile2 [open rmsd_r2_$no.csv w]
            #set outfile3 [open rmsd_a1_$no.csv w]
            #set outfile4 [open rmsd_a2_$no.csv w]
            
            for {set frame $minFrame} {$frame < $maxFrame} {incr frame} {
                
                $comparison1 frame $frame
                $comparison2 frame $frame
                #$comparison3 frame $frame
                #$comparison4 frame $frame

                set rmsd1 [measure rmsd $comparison1 $reference1]
                set rmsd2 [measure rmsd $comparison2 $reference2]
                #set rmsd3 [measure rmsd $comparison3 $reference3]
                #set rmsd4 [measure rmsd $comparison4 $reference4]
                
                set adframe [expr $frame - 499]
                puts $outfile1 "$adframe, $rmsd1"
                puts $outfile2 "$adframe, $rmsd2"
                #puts $outfile3 "$adframe, $rmsd3"
                #puts $outfile4 "$adframe, $rmsd4"
                
            }
            
            $reference1 delete
            $reference2 delete
            #$reference3 delete
            #$reference4 delete
            $comparison1 delete
            $comparison2 delete
            #$comparison3 delete
            #$comparison4 delete
            
            close $outfile1
            close $outfile2
            #close $outfile3
            #close $outfile4
            
        }
    }
}