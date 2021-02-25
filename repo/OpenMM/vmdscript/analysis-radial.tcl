proc rmsd_radial_resting {start R1 xori yori num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set Rsq1 [expr $R1*$R1]
             
             set sel [atomselect top "name RM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open restingrmsd$no.csv w]
             
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

proc radial_density_resting {start R1 R2 xori yori num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set Rsq1 [expr $R1*$R1]
            set Rsq2 [expr $R2*$R2]

            set outfile [open resting_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                set sel1 [atomselect top "name RM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) < $Rsq1" frame $f]
                set sel2 [atomselect top "name RM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq1 and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) < $Rsq2 " frame $f]
                set sel3 [atomselect top "name RM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq2 " frame $f]
                                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                
                puts $outfile "$f, $num1, $num2, $num3"
                
                unset num1
                unset num2
                unset num3
                
                unset sel1
                unset sel2
                unset sel3
            }
            close $outfile
        }
        
    }
}

proc rmsd_radial_activated {start R1 xori yori num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series} {incr k} {
         for {set i 0} {$i < $repeat} {incr i} {
             set no [expr {$k + $num_sim_series*$i }]
             mol load pdb test$no.pdb dcd test$no.dcd
             set Rsq1 [expr $R1*$R1]
             
             set sel [atomselect top "name AM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq1" frame $start]
             set index1 [ lsort -integer -index 0 -increasing -unique [$sel get {index}]]
             set outfile [open activatedrmsd$no.csv w]
             
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
 

proc radial_density_activated {start R1 R2 xori yori num_sim_series repeat} {
     for {set k 1 } {$k <= $num_sim_series } {incr k} {
         for {set i 0 } {$i < $repeat } {incr i} {
            set no [expr {$k + $num_sim_series*$i }]
            mol load pdb test$no.pdb dcd test$no.dcd
            set Rsq1 [expr $R1*$R1]
            set Rsq2 [expr $R2*$R2]

            set outfile [open activated_density$no.csv w]
             
            for {set f $start } {$f < [molinfo top get numframes]} {incr f} {
                set sel1 [atomselect top "name AM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) < $Rsq1" frame $f]
                set sel2 [atomselect top "name AM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq1 and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) < $Rsq2 " frame $f]
                set sel3 [atomselect top "name AM and (x-$xori)*(x-$xori)+(y-$yori)*(y-$yori) > $Rsq2 " frame $f]
                                
                set num1 [$sel1 num]
                set num2 [$sel2 num]
                set num3 [$sel3 num]
                
                puts $outfile "$f, $num1, $num2, $num3"
                
                unset num1
                unset num2
                unset num3
                
                unset sel1
                unset sel2
                unset sel3
            }
            close $outfile
        }
        
    }
}