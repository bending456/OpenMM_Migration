puts "this script is for measuring distance that each particle travels within a given simulation time"
proc stats { minFrame maxFrame selection } {
    set sel0 [atomselect top $selection]
    set index_list [lsort -integer -index 0 -increasing -unique [$sel0 get {index}]]
    set no_of_cells 0
    set sum_x_dist 0
    set sum_y_dist 0

    foreach particle $index_list {
        set no_of_cells [expr $no_of_cells + 1]
        set sel1 [atomselect top "index $particle"]
        $sel1 frame $minFrame
        
        set coordmin [measure center $sel1] 
        set xcoordmin [lindex $coordmin 0]
        set ycoordmin [lindex $coordmin 1]
               
        $sel1 frame $maxFrame
        set coordmax [measure center $sel1]
        set xcoordmax [lindex $coordmax 0]
        set ycoordmax [lindex $coordmax 1]

        set x_dist [expr abs($xcoordmax - $xcoordmin)]
        set y_dist [expr abs($ycoordmax - $ycoordmin)]

        set sum_x_dist [expr $sum_x_dist + $x_dist]
        set sum_y_dist [expr $sum_y_dist + $y_dist]

    }
    
    set x_avg_dist [expr $sum_x_dist/$no_of_cells]
    set y_avg_dist [expr $sum_y_dist/$no_of_cells]

    set x_diff 0
    set y_diff 0

    foreach particle $index_list {
        set sel2 [atomselect top "index $particle"]
        $sel2 frame $minFrame

        set coordmin [measure center $sel2] 
        set xcoordmin [lindex $coordmin 0]
        set ycoordmin [lindex $coordmin 1]
               
        $sel2 frame $maxFrame
        set coordmax [measure center $sel2]
        set xcoordmax [lindex $coordmax 0]
        set ycoordmax [lindex $coordmax 1]

        set x_dist [expr abs($xcoordmax - $xcoordmin)]
        set y_dist [expr abs($ycoordmax - $ycoordmin)]

        set x_diff [expr ($x_diff + pow(($x_dist-$x_avg_dist),2))]
        set y_diff [expr ($y_diff + pow(($y_dist-$y_avg_dist),2))]
        
    }

    set x_std [expr sqrt($x_diff/$no_of_cells)]
    set y_std [expr sqrt($y_diff/$no_of_cells)]
    set x_serr [expr $x_std/sqrt($no_of_cells)] 
    set y_serr [expr $y_std/sqrt($no_of_cells)]

    puts "$no_of_cells"

    puts "Avg Dist. in X-dir $x_avg_dist"
    puts "with standard deviation of $x_std"
    puts "with standard error of $x_serr"
    puts "Avg Dist. in Y-dir $y_avg_dist"
    puts "with standard deviation of $y_std"
    puts "with standard error of $y_serr"

}

proc multiple { minFrame maxFrame sel1 sel2 } {
    stats $minFrame $maxFrame $sel1 
    stats $minFrame $maxFrame $sel2
}