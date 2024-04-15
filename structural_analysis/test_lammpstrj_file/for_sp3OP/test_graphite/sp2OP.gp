set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "sp2OP.pdf"

#set xrange [0:70]
set yrange [0:0.5]

set xlabel "X (Angstrom)"
set ylabel "Hexagonal (sp2) order parameter"

p 'sp2OP_profile_graphite_288.lammpstrj.dat' u 1:2 w p pt 5 ps 0.25 lc rgb "#0000FF" title ""

