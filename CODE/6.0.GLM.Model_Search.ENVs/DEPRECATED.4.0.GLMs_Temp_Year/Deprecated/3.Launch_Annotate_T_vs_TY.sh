#!/bin/bash
#SBATCH -J annotate_glm # A single job name for the array
#SBATCH --ntasks-per-node=4 
#SBATCH -N 1 # on one node
#SBATCH -t 3:00:00 
#SBATCH --mem 16G
#SBATCH -o /scratch/yey2sn/Overwintering_ms/4.GML_plots/enrrich_slurm/annot.%A_%a.out # Standard output
#SBATCH -e /scratch/yey2sn/Overwintering_ms/4.GML_plots/enrrich_slurm/annot.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account=jcbnunez
#SBATCH --array=1-1232

echo "begin loop"

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1

guide=./3.annotation_guide_file.txt

start=$( cat $guide | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
stop=$( cat $guide  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )

echo $start $stop

Rscript --vanilla 3.Annotate_T_vs_TY.r $start $stop