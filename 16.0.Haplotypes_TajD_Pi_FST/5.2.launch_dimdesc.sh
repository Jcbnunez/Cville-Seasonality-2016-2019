#!/usr/bin/env bash
#
#SBATCH -J run_dimdesc # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 6 hours
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1

##
Rscript --vanilla 5.1.run.dim.desc.r

