#!/bin/bash -l

set -o pipefail

#$ -S /bin/bash
#$ -l h_rt=12:00:00
#$ -l mem=20G
#$ -l tmpfs=50G
#$ -cwd


number=$SGE_TASK_ID
paramfile=${PWD}/run_${RAND_ID}.txt
 
run=`sed -n ${number}p $paramfile`
$run

