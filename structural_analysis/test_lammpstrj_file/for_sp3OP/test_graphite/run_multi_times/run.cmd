#!/bin/bash

code='/usr/workspace/wsa/lyu1/utils/structural_analysis/OP'

for num in 1 2 3 4 5 6 7 8 9 10
do

  $code graphite_288.lammpstrj 1.7 12 12 1.0 1 > test.out

  mv test.out test_$num.out

  mv OP_profile_graphite_288.lammpstrj.dat OP_profile_graphite_288.lammpstrj_$num.dat
	
done
