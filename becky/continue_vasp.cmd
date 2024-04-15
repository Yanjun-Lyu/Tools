#!/bin/bash
#MSUB -N dntf_383
#MSUB -l nodes=10:ppn=36
#MSUB -l walltime=24:00:00
#MSUB -q pbatch
#MSUB -A carbnuc
#MSUB -V 
#MSUB -m abe
#MSUB -o stdoutmsg

module load mkl

last_run=0

# Find the index used to set file names

if [ -e ${last_run}.INCAR ] ; then

	# This is a continuation job... figure out which
	
	last_run=`ls *.INCAR | tail -n 1` # Something like "4.INCAR"
	last_run=${last_run%*.INCAR}
	
	# Increment
	
	let last_run=last_run+1	
fi

for i in INCAR POSCAR OUTCAR OSZICAR CONTCAR stdoutmsg
do
	cp $i ${last_run}.${i}
done

cp CONTCAR POSCAR
	

srun -N 10 -n 360 /usr/gapps/emc-vasp/vasp.5.4.1/build/std/vasp > ${TAG}.out  
