#!/bin/bash
#SBATCH -t 6-20:00:00
#SBATCH -N 1
#SBATCH --mem=9G
#SBATCH -p standard
#SBATCH -c 16
#SBATCH -A berglandlab_standard

#WIP

module load gcc/9.2.0
module load openmpi/3.1.6
module load R/4.2.1
baypass="/home/nzx3cc/baypass_2.4/sources/g_baypass"
$baypass -gfile /scratch/nzx3cc/nzx3cc/rawdata_baypass/wholegenome.genobaypass -poolsizefile /scratch/nzx3cc/nzx3cc/rawdata_baypass/wholegenome.poolsize -outprefix "allchr" -nthreads 16 -omegafile /scratch/nzx3cc/nzx3cc/rawdata_baypass/allthinned_mat_omega.out -efile /scratch/nzx3cc/nzx3cc/env_analysis/envdata.txt -contrastfile /scratch/nzx3cc/nzx3cc/env_analysis/contrasttable.txt

