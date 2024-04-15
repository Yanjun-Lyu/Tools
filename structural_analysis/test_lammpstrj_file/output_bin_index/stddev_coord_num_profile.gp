set terminal pdf  size 2.5,2.0 enhanced color font "Times,12" background rgb 'white'
set encoding iso_8859_1 
set output "stddev_coord_num_profile.pdf"

set xrange [0:70]
#set yrange [3:4.2]

set xlabel "X (\AA)"
set ylabel "Standard deviation: coordination number"

z(x) = (x - 0.5) * 0.5

p 'stddev_coord_num_profile.dat'          u (z($1)):($2) w p pt 5 ps 0.25 lc rgb "#DD0000FF" title "", \
  'stddev_coord_num_profile_tanh_fit.dat' u (z($1)):($2) w l lc rgb "#FF0000" title ""

