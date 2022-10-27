#!/usr/bin/env bash
#
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=85G
#SBATCH --time=36:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

### load modules
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 
module load samtools htslib bcftools/1.9 gparallel/20170822
module load tabix

###########################################################################

shapeit_out=/scratch/yey2sn/phasing_droso_data/Whatshap/shapeit_out

PIPELINE=CM_pops

###########################################################################


## Part 1. Make a list of the bcf files
ls -d $shapeit_out/$PIPELINE.*.shapeit.bcf \
 > $shapeit_out/bcfs.$PIPELINE.list

## Part 2. Index the bcf tools
for f in $shapeit_out/$PIPELINE.*.shapeit.bcf; do
  echo "Index file -> $f"
  bcftools index \
  -f \
  ${f}
done

# Part 3. Concatenate all bcfs
bcftools \
concat \
-f $shapeit_out/bcfs.$PIPELINE.list \
-O v \
-l \
-o $shapeit_out/$PIPELINE.AllChrs.Whatshap.shapeit.vcf

# Part 4. bgzip and tabix
bgzip $shapeit_out/$PIPELINE.AllChrs.Whatshap.shapeit.vcf
tabix $shapeit_out/$PIPELINE.AllChrs.Whatshap.shapeit.vcf.gz

 echo "done" $(date)
