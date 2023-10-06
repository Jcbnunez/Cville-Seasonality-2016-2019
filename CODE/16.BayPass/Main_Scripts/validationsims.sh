#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -A berglandlab
#SBATCH -c 16
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem 9G
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p podcontrol.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt2
echo $opt1
$baypass -gfile G.pod_${opt1} -outprefix tempmaxpod_${opt1}_${opt2} -nthreads 16 -omegafile /scratch/nzx3cc/nzx3cc/rawdata_baypass/allthinned_mat_omega.out -efile /scratch/nzx3cc/nzx3cc/env_analysis/tmax.txt
