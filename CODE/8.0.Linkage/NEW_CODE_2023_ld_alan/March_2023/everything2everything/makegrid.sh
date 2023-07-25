#!/bin/bash
#
#SBATCH -J runGrid # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:60:00 ### 30 min
#SBATCH --mem 100G
#SBATCH -o /scratch/aob2x/everythingLD/logs/runGrid.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/everythingLD/logs/runGrid.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard



### sbatch --array=1-1000 ~/Overwintering_18_19/ld_alan/everything2everything/makegrid.sh
### sbatch --array=15,18,31,36,39,45 ~/Overwintering_18_19/ld_alan/everything2everything/makegrid.sh
### sbatch --array=15,63,64,68,83,92,101,102,105,125,126,150,151,155,156,158,159,160,161,169,173,184,199,200,201,209,212,213,246,255,258,268,292,308,353,486,487,498 ~/Overwintering_18_19/ld_alan/everything2everything/makegrid.sh
### sacct -j 49921682
### cat /scratch/aob2x/everythingLD/logs/runGrid.49889342*.out
### cat /scratch/aob2x/everythingLD/logs/runGrid.49889342*.err

#sacct -j 49906492 | grep "TIME" | cut -f2 -d'_' | cut -f1 -d' ' | tr '\n' ','

module load intel/18.0 intelmpi/18.0 R/3.6.3

#SLURM_ARRAY_TASK_ID=1
date

cd ~
Rscript Overwintering_18_19/ld_alan/everything2everything/makegrid.R ${SLURM_ARRAY_TASK_ID} 1000

date
