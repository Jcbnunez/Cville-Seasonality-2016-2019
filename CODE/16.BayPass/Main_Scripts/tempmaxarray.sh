#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 9G
#SBATCH -t 4-20:30:00
#SBATCH -p standard
#SBATCH -A berglandlab
a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p tempmaxcontrol.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt2
echo $opt1
echo $opt3

baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
  $baypass -gfile /scratch/nzx3cc/nzx3cc/rawdata_baypass/"${opt2}".genobaypass -poolsizefile /scratch/nzx3cc/nzx3cc/rawdata_baypass/${opt2}.poolsize -outprefix "tempmax_${opt2}_${opt3}" -nthreads 16  -omegafile /scratch/nzx3cc/nzx3cc/rawdata_baypass/"${opt1}"thinned_mat_omega.out -efile /scratch/nzx3cc/nzx3cc/env_analysis/tmax.txt 

 
