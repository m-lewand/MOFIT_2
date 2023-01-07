set terminal gif animate delay 20 size 1000,1000
set output 'Eigenstate.gif'
stats 'Eigenstates_time.dat' nooutput
set xrange [-40:40]
set yrange [-40:40]
set view map
do for [i=1:int(STATS_blocks)] {
    sp 'Eigenstates_time.dat' u 1:2:3 index (i-1) palette z pt 7 ps 2
}
