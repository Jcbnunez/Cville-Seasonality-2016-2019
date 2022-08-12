#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40G
#SBATCH --time=32:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

# This script will remove repetitive regions and low complexity regions form the VCF

module load gcc/9.2.0
module load bedtools/2.29.2
module load tabix
module load bcftools

#Declare variables
IN_GZVCF=/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.vcf.gz
repetitive_regions=/project/berglandlab/Dmel_genomic_resources/Repeats/Dmel6_combined_repeats.sort.merge.bed
out_name=dest.all.PoolSNP.001.50.10Mar2021.ann.noRep

#Remove options
bedtools intersect \
-v -header \
-a $IN_GZVCF \
-b <(sort -k1,1 -k2,2n $repetitive_regions) \
> $out_name.vcf

#bgzip and tabix
bgzip $out_name.vcf
tabix $out_name.vcf.gz