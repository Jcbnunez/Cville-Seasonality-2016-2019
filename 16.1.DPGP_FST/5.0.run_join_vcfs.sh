#!/usr/bin/env bash
#
#SBATCH -J merge.DPGP # A single job name for the array
#SBATCH --ntasks-per-node=2 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 15:00:00 ### 
#SBATCH --mem=120G 
#SBATCH -o slurmOut/merge.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/merge.ag%A_%a.err # Standard error
#SBATCH --partition=largemem
#SBATCH --account berglandlab

#files.to.merge.decompressed.list --> remove bz2 to make .vcf 
# also add the hard address --> /scratch/yey2sn/Supergene_paper/2.DPGP/individual.vcfs/ZI85_sites.vcf_fold/ZI85_sites_full.vcf

module load bcftools
 
bcftools merge \
--missing-to-ref \
vcfs_diploidized_header/*.vcf.gz \
> DPGP3.2L.merged.vcf

bgzip DPGP3.2L.merged.vcf
tabix DPGP3.2L.merged.vcf.gz

echo done