set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "coord_num_profile.pdf"

set xrange [0:70]
set yrange [3:4.2]

set xlabel "X (\AA)"
set ylabel "Average coordination number"

z(x) = (x - 0.5) * 0.5

p 'coord_num_profile_traj_30GPa_5500K_last_100_steps.lammpstrj.dat' u (z($1)):($2) w p pt 5 ps 0.25 lc rgb "#DD0000FF" title ""

