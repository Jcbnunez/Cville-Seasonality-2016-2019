#!/usr/bin/env bash
#
#SBATCH -J std.calc_ihs # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 
#SBATCH --mem 80G
#SBATCH -p standard
#SBATCH --account jcbnunez


module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1

##
Rscript --vanilla 3.calculate_ihs_STD.r
