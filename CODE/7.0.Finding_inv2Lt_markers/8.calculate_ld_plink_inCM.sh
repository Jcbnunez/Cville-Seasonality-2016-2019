#!/usr/bin/env bash
#
#SBATCH -J ld_plink_CM # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 4 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut/ld_plink_cm.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/ld_plink_cm.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#load module
module load plink
module load bcftools

#make folder
out_folder=$1
mkdir $out_folder

#load cpu info
cpu=$SLURM_CPUS_PER_TASK

#load input data
vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.wSNPids.vcf.gz
snp_guide_file=/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/Inv2L_markers_to_use_CM.txt

#declare array variables
snp=$(cat $snp_guide_file | sed '1d'  | cut -f 7,11 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 2 )
rank=$(cat $snp_guide_file | sed '1d'  | cut -f 7,11 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 1 )

echo $snp
echo $rank

bcftools query -l $vcf

#run plink
plink --vcf $vcf \
--allow-extra-chr \
--r2 \
--double-id \
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




