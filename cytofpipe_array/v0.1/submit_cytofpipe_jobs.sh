#!/bin/bash -l

#$ -S /bin/bash
#$ -l h_rt=12:00:00
#$ -l mem=20G
#$ -l tmpfs=50G
#$ -cwd


module unload compilers
module unload mpi

#– this was the old r/recommended
module load gcc-libs            
module load compilers/gnu/4.9.2 
module load openblas/0.2.14/gnu-4.9.2
module load mpi/openmpi/1.10.1/gnu-4.9.2
module load java/1.8.0_92
module load fftw/3.3.4/gnu-4.9.2
module load ghostscript/9.19/gnu-4.9.2
module load texinfo/5.2/gnu-4.9.2
module load texlive/2015
module load gsl/1.16/gnu-4.9.2
module load hdf/5-1.8.15/gnu-4.9.2
module load netcdf/4.3.3.1/gnu-4.9.2
module load udunits/2.2.20/gnu-4.9.2
module load jags/4.2.0/gnu.4.9.2-openblas
module load root/5.34.30/gnu-4.9.2
module load glpk/4.60/gnu-4.9.2
module load perl/5.22.0
module load libtool/2.4.6
module load graphicsmagick/1.3.21
module load python2/recommended
module load gdal/2.0.0
module load gmt/recommended
module load proj.4/4.9.1
module load geos/3.5.0/gnu-4.9.2
module load protobuf/12-2017/gnu-4.9.2
module load jq/1.5/gnu-4.9.2

module load r/3.4.2-openblas/gnu-4.9.2


number=$SGE_TASK_ID
paramfile=${PWD}/run_${RAND_ID}.txt
 
run=`sed -n ${number}p $paramfile`
$run

