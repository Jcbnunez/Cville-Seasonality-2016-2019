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
module load gcc/9.2.0 shapeit4
module load vcftools
module load tabix

# import guide files
# All VCF files must be bgzipped and tabix'ed. Including the DGRP.

natural_pops=/scratch/yey2sn/phasing_droso_data/Whatshap/natural_pops_sampleid.txt
intervals=/scratch/yey2sn/phasing_droso_data/Whatshap/Intervals_for_phasing.txt
MERGED_WHATSHAP=/scratch/yey2sn/phasing_droso_data/Whatshap/MERGED_WHATSHAP_noSRs
DGRP=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.vcf.gz
PIPELINE=CM_pops

WD=/scratch/yey2sn/phasing_droso_data/Whatshap/
## make sure we are in the WD
cd $WD

###########################################################################

chr=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`
echo $chr

###########################################################################
# Generate folders

if [ -d "shapeit_out" ]
then
	echo "Working mapping_stats folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir shapeit_out
	date
fi

###########################################################################

# Part 1: Phase only individuals from natural populations

### Filter

 vcftools \
	--gzvcf $MERGED_WHATSHAP/${chr}.whatshapp.noSR.vcf.gz \
	--keep $natural_pops \
	--recode \
	--recode-INFO-all \
	--out ./shapeit_out/$PIPELINE.${chr}

### Tabix and indes

bgzip ./shapeit_out/$PIPELINE.${chr}.recode.vcf
tabix ./shapeit_out/$PIPELINE.${chr}.recode.vcf.gz

### run shapeit

shapeit4 \
  --input ./shapeit_out/$PIPELINE.${chr}.recode.vcf.gz \
  --region ${chr} \
  --use-PS 0.0001 \
  --thread $CPU \
  --seed 123456 \
  --reference $DGRP \
  --sequencing \
  --log ./shapeit_out/${chr}.daphnid.log \
  --ibd2-length 5 \
  --ibd2-output ./shapeit_out/${chr}.onPerSC.daphnid.IBD2blacklist.txt.gz \
  --output ./shapeit_out/$PIPELINE.${chr}.shapeit.bcf
  
 echo "done" $(date)
