set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "OP.pdf"

set xrange [-5:65]
set yrange [0:0.3]

set xlabel "X (Angstrom)"
set ylabel "Order parameters"

p 'OP_profile_traj_100GPa_6450GPa_last_frame.lammpstrj.dat' u 1:2 w lp pt 7 ps 0.1 lc rgb "#FF0000" title "q_{tet}", \
  'OP_profile_traj_100GPa_6450GPa_last_frame.lammpstrj.dat' u 1:3 w lp pt 7 ps 0.1 lc rgb "#0000FF" title "q_{hex}"

