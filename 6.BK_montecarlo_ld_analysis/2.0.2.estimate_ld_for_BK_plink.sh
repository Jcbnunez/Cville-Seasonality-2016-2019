#!/bin/bash
#SBATCH -J ld_plink_bk # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 2 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut/ld_plink_dgrp.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/ld_plink_dgrp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account jcbnunez

#load module
module load plink

mkdir bk_ld_output

#Input 1 = "type of SNP as the focus snps"
#Input 2 = "the SNP file to iterate over"
#Input 3 = "the focal snp"

snp_guide_file=/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/guide_files_bk_ld.txt

#load information
focal_type=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 2 )  #----------> "inv_focus"
#load list of SNPs to iterate 
iterator_snps_list=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 3 ) # --> /scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/temperature_snps_ids_top100.txt
#set a fix marker SNP
focal_snp=$(cat $snp_guide_file  | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 1 )  # --->2L_12899672_SNP

echo $focal_type
echo $iterator_snps_list
echo $focal_snp

# fix inputs
#load input data
vcf=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/ForPlinkMarkersOnly.CM.2L.recode.vcf.gz

#loop over temperature outlier snps

while read dependent_snp
#dependent_snp=2L_748285_SNP
do # <---------- do here

#run plink
plink --vcf $vcf \
--allow-extra-chr \
--double-id \
--ld $focal_snp $dependent_snp \
--out tmp.$focal_snp.$dependent_snp

#clean up
rm tmp.$focal_snp.$dependent_snp.nosex

#slice out the valuable information
cat tmp.$focal_snp.$dependent_snp.log | \
sed '1,33d' | \
head -n 4 | \
sed -E 's/\ +/\ /g' | \
sed -E 's/\ /\t/g' | \
sed -n '1p;2p' > tmp.$focal_snp.$dependent_snp.flt.txt

rm tmp.$focal_snp.$dependent_snp.log

if [ -s tmp.$focal_snp.$dependent_snp.flt.txt ]
then

paste -d "\t" tmp.$focal_snp.$dependent_snp.flt.txt <(printf %"s\n" "present" "absent") \
> tmp.$focal_snp.$dependent_snp.flt.hdr.txt

rm tmp.$focal_snp.$dependent_snp.flt.txt

sed -i "s/^/${focal_snp}vs${dependent_snp}/" tmp.$focal_snp.$dependent_snp.flt.hdr.txt
sed -i "s/^/$focal_type\t/" tmp.$focal_snp.$dependent_snp.flt.hdr.txt

cat tmp.$focal_snp.$dependent_snp.flt.hdr.txt >> bk_ld_output/$focal_snp.marker.BK.txt
rm tmp.$focal_snp.$dependent_snp.flt.hdr.txt

else
	 echo $dependent_snp "does not exist"
	 rm tmp.$focal_snp.$dependent_snp.flt.txt
fi

done < $iterator_snps_list # <---------- done
