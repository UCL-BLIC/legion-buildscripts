#!/bin/bash -l

#$ -S /bin/bash
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -cwd


module unload compilers
module unload mpi
module load r/recommended


number=$SGE_TASK_ID
paramfile=${PWD}/run_${RAND_ID}.txt
 
run=`sed -n ${number}p $paramfile`
$run

