#!/usr/bin/env bash
#
#SBATCH -J remake.BPGP # A single job name for the array
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

module load tabix 

loc_i=$( cat 2.files.to.merge.decompressed.list | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{print $1}' )
echo $loc_i

bgzip $loc_i
tabix $loc_i.gz

echo "done"
date