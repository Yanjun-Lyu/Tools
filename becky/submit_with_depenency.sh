#!/bin/bash

curr[0]=3225452 # /p/lscratchh/rlindsey/DNTF_CALCULATIONS/1.80gcc_383K
curr[1]=3225455 # /p/lscratchh/rlindsey/DNTF_CALCULATIONS/2.50gcc_9000K

idx=0

for i in 1.80gcc_383K  2.50gcc_9000K
do
	cd $i
	
	jobno=0
	
	echo "On state point $i:"

	for j in {1..20} # Number of resumbissions 
	do

		if [ $j -eq 1 ] ; then
			jobno=${curr[$idx]}
		fi
		
		echo "	Setting jobno: $jobno "
		
		awk -v job="#MSUB -l depend=${jobno}" '/carbnuc/{print;getline;print job}!/carbnuc/{print}' continue_vasp.cmd > continue_vasp-dependency_${jobno}.cmd
		
		echo "		Saved to file continue_vasp-dependency_${jobno}.cmd"
		
		tmp_jobno=`msub continue_vasp-dependency_${jobno}.cmd`
		tmp_jobno=`echo $tmp_jobno | awk '{print $1}'`

		jobno=$tmp_jobno
		
		echo "			...submitted job $jobno"

	done
	
	cd ..
	
	let idx=idx+1
done

		
