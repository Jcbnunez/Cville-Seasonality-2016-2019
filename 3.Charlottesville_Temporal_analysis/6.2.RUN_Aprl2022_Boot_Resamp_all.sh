#!/usr/bin/env bash
#
#SBATCH -J resample # A single job name for the array
#SBATCH -c 3
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 #<= this may depend on your resources
#SBATCH --mem 12G #<= this may depend on your resources
#SBATCH -o ./slurmOut/resample.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/resample.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH -A jcbnunez
#SBATCH --array=1-7500

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj


guide=guide_file_2022_boot_resamp.28Cov.varCov.txt

jobid=$( cat $guide  |  awk '{ print $1 }' |  sed "${SLURM_ARRAY_TASK_ID}q;d" )
snpnum=$( cat $guide |  awk '{ print $2 }' |  sed "${SLURM_ARRAY_TASK_ID}q;d" )
cov=$( cat $guide    |  awk '{ print $3 }' |  sed "${SLURM_ARRAY_TASK_ID}q;d" )
 
echo $jobid
echo $snpnum
echo $cov

Rscript \
--vanilla \
6.0.Aprl2022_Boot_Resamp_all.R \
$jobid $snpnum $cov