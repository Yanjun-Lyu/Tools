set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "sp3OP.pdf"

#set xrange [0:70]
#set yrange [3:4.2]

set xlabel "X (Angstrom)"
set ylabel "Tetrahedral (sp3) order parameter"

p 'sp3OP_profile_diamond_64.lammpstrj.dat' u 1:2 w p pt 5 ps 0.25 lc rgb "#0000FF" title ""

