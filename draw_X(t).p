	set terminal png size 1000,1000
	set output 'X(t).png'
	set lmargin at screen 0.2
	set rmargin at screen 0.8
	set xlabel 't [atomic units]' font ",25"
	set ylabel 'X [nm]' font ",25"
	plot 'X.dat' w l lw 2
