set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "OP.pdf"

set key at 8, 0.5

#set xrange [0:70]
set yrange [-0.1:1.0]

set xlabel "Z (Angstrom)"
set ylabel "Order parameters"

p 'z_OP_profile_graphite_288.lammpstrj.dat' u 1:2 w lp pt 7 ps 0.1 lc rgb "#FF0000" title "q_{tet}", \
  'z_OP_profile_graphite_288.lammpstrj.dat' u 1:3 w lp pt 7 ps 0.1 lc rgb "#0000FF" title "q_{hex}"

