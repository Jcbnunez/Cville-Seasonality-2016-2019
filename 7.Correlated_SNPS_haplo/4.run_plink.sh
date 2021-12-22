#!/usr/bin/env bash
#
#SBATCH -J ld_plink # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 4 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut_ld/ld_plink_dgrp.%A_%a.out # Standard output
#SBATCH -e ./slurmOut_ld/ld_plink_dgrp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account jcbnunez
#SBATCH --array=1-3446%50

### 3446

#load module
module load plink

#make folder
out_folder=ld_plink
mkdir $out_folder

#load cpu info
cpu=$SLURM_CPUS_PER_TASK

#load input data
vcf=vcf_for_ld_analysis.2l.vcf.recode.vcf.gz
snp_guide_file=SNPS_to_select.txt

#declare array variables
snp=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d"  )

echo $snp

#run plink
plink --vcf $vcf \
--allow-extra-chr \
--double-id \
--r2 \
--ld-snp $snp \
--ld-window 99999999 \
--ld-window-kb 99999 \
--threads $cpu \
--ld-window-r2 0.0 \
--out $out_folder/$snp.ld

#clean up
rm $out_folder/$snp.ld.nosex
rm $out_folder/$snp.ld.log

