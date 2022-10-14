#!/usr/bin/env bash
#
#SBATCH -J fix.BPGP # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 
#SBATCH --mem 20G
#SBATCH -o slurmOut/glm.fsts.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/glm.fsts.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-201

#files.to.merge.decompressed.list --> remove bz2 to make .vcf 
# also add the hard address --> /scratch/yey2sn/Supergene_paper/2.DPGP/individual.vcfs/ZI85_sites.vcf_fold/ZI85_sites_full.vcf

mkdir fixedGTs

module load tabix 
module load vcftools

loc_i=$( cat 2.files.to.merge.decompressed.list | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{print $1}' )
echo $loc_i
name_i=$( cat 2.files.to.merge.decompressed.list | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{print $2}' )
echo $name_i


vcftools --gzvcf $loc_i.gz \
--remove-indels \
--min-alleles 1 \
--max-missing  1 \
--recode --recode-INFO-all \
--out fixedGTs/$name_i

sed 's|PL\t0:|PL\t0/0:|g' fixedGTs/$name_i.recode.vcf > fixedGTs/$name_i.recode.fixGT.vcf

bgzip fixedGTs/$name_i.recode.fixGT.vcf

tabix fixedGTs/$name_i.recode.fixGT.vcf.gz


echo "done"
date