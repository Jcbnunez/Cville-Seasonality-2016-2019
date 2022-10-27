#!/usr/bin/env bash
#
#SBATCH -J ld_plink_dgrp # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 4 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut/ld_plink_dgrp.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/ld_plink_dgrp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#load module
module load plink

#make folder
out_folder=$1
mkdir $out_folder

#load cpu info
cpu=$SLURM_CPUS_PER_TASK

#load input data
vcf=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz
keep_ids=DGRP_lines_used_for_2lt_demarcation_plink.txt
snp_guide_file=inv2L_informative_markers_Dm3.txt

#declare array variables
snp=$(cat $snp_guide_file | sed '1d'  | cut -f 3,7 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 1 )
rank=$(cat $snp_guide_file | sed '1d'  | cut -f 3,7 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 2 )

echo $snp
echo $rank


#run plink
plink --vcf $vcf \
--allow-extra-chr \
--keep $keep_ids \
--r2 \
--ld-snp $snp \
--ld-window 99999999 \
--ld-window-kb 99999 \
--threads $cpu \
--ld-window-r2 0.1 \
--out $out_folder/$snp.$rank

#add rank information
sed -i "s/^/$snp|$rank\ /" $out_folder/$snp.$rank.ld

#clean up
rm $out_folder/$snp.$rank.nosex
rm $out_folder/$snp.$rank.log

