set terminal gif animate delay 10
set output 'foobar.gif'
stats 'Eigenstates_time.dat' nooutput
set xrange [-50:50]
set yrange [-50:50]
set view map
do for [i=1:int(STATS_blocks)] {
    sp 'Eigenstates_time.dat' u 1:2:3 index (i-1) palette z pt 7 ps 2
}
