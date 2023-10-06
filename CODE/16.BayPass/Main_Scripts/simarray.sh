#!/bin/bash
#SBATCH -t 6-20:00:00
#SBATCH -A berglandlab
#SBATCH -c 16
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem 9G 
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
Rscript /scratch/nzx3cc/nzx3cc/scripts/simscript.R --args "pod_${SLURM_ARRAY_TASK_ID}" all 577000
$baypass -gfile G.pod_${SLURM_ARRAY_TASK_ID} -outprefix anapod_${SLURM_ARRAY_TASK_ID} -nthreads 8 -omegafile /scratch/nzx3cc/nzx3cc/rawdata_baypass/allthinned_mat_omega.out -efile /scratch/nzx3cc/nzx3cc/env_analysis/envdata.txt

