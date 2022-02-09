#!/bin/bash
#SBATCH -J ld_plink_bk # A single job name for the array
#SBATCH -c 1 
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 2 hours
#SBATCH --mem 4G
#SBATCH -o ./slurmOut/ld_plink_dgrp.%A_%a.out # Standard output
#SBATCH -e ./slurmOut/ld_plink_dgrp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account jcbnunez

#load module
module load plink

mkdir bk_ld_output

### import guide file
snp_guide_file=$1

#load focal SNP
focal_snp=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 1 )  
#load affected SNP
#affect_snp=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 3 )  
#load comparison type
#comp_type=$(cat $snp_guide_file | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f 5 )  

echo $focal_snp
#echo $affect_snp
#echo $comp_type

# fix inputs
#load input data
vcf=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/ForPlinkMarkersOnly.CM.2L.recode.vcf.gz

#loop over temperature outlier snps
while read affect_snp
do
#run plink
plink --vcf $vcf \
--allow-extra-chr \
--double-id \
--ld $focal_snp $affect_snp \
--out tmp.$focal_snp.$affect_snp

#clean up
rm tmp.$focal_snp.$affect_snp.nosex

#slice out the valuable information
cat tmp.$focal_snp.$affect_snp.log | \
sed '1,33d' | \
head -n 4 | \
sed -E 's/\ +/\ /g' | \
sed -E 's/\ /\t/g' | \
sed -n '1p;2p' > tmp.$focal_snp.$affect_snp.flt.txt

rm tmp.$focal_snp.$affect_snp.log

if [ -s tmp.$focal_snp.$affect_snp.flt.txt ]
then

paste -d "\t" tmp.$focal_snp.$affect_snp.flt.txt <(printf %"s\n" "present" "absent") \
> tmp.$focal_snp.$affect_snp.flt.hdr.txt

rm tmp.$focal_snp.$affect_snp.flt.txt

sed -i "s/^/${focal_snp}|${affect_snp}/" tmp.$focal_snp.$affect_snp.flt.hdr.txt

cat tmp.$focal_snp.$affect_snp.flt.hdr.txt >> bk_ld_output/$focal_snp.Output.marker.BK.txt
rm tmp.$focal_snp.$affect_snp.flt.hdr.txt

else
	 echo $dependent_snp "does not exist"
	 rm tmp.$focal_snp.$affect_snp.flt.txt
fi

done <$snp_guide_file

date