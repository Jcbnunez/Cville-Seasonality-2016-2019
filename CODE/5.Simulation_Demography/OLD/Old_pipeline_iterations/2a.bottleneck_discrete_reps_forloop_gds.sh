#!/usr/bin/env bash
#
#SBATCH -J slim-bottleneck # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-06:00 #  hours
#SBATCH --mem 16G
#SBATCH -o /scratch/csm6hg/slim_bottleneck/err/bottle.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/slim_bottleneck/err/bottle.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

SLURM_JOB_ID=1
SLURM_ARRAY_TASK_ID=503

# Working directory
wd="/scratch/csm6hg/slim_bottleneck"

# Parameter file
paramFile=${wd}/model_paramList6

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
nMax=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
nMin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )

# Output directory
output="/project/berglandlab/connor/BACKUP_scratch/slim_bottleneck/data"

# For loop replicate
for replicate in {1..100}; do

  # Progress
  echo 'Replicate:' $replicate

  # Random seed
  seed=$(($RANDOM))

  # RAM disk set up
  [ ! -d ${wd}/tmp/ ] && mkdir ${wd}/tmp/
  [ ! -d ${wd}/tmp/${SLURM_JOB_ID} ] && mkdir ${wd}/tmp/${SLURM_JOB_ID}
  [ ! -d ${wd}/tmp/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir ${wd}/tmp/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
  [ ! -d ${wd}/tmp/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${seed} ] && mkdir ${wd}/tmp/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${seed}

  # Temporary directory for VCFs and metadata
  tmpdir=${wd}/tmp/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${seed}

  # Load modules
  module load samtools/1.10 parallel
  module load htslib/1.9
  module load gcc/9.2.0 slim

  # Run SLiM
  slim \
  -d slurmIDi=$slurmID \
  -d StartPopi=1500 \
  -d nMaxi=$nMax \
  -d nMini=$nMin \
  -d Repi=$Rep \
  -d nSampi=$nSamp \
  -d Geni=$Gen \
  -d simIDi=$seed \
  -d SJi=$SLURM_JOB_ID \
  -d STi=$SLURM_ARRAY_TASK_ID \
  -d replicatei=${replicate} \
  ${wd}/1a.bottleneck_discrete.slim

  # Creates list for all VCFs
  ls -d ${tmpdir}/*vcf > \
  ${tmpdir}/vcf.list

  # gzip and tabix vcf files
  find ${tmpdir} -maxdepth 1 -type f -exec bgzip {} \;
  find ${tmpdir} -maxdepth 1 -type f -exec tabix {} \;

  # Run through Rscript
  echo "Running through summary scripts"

  # Load R module
  module load goolf/7.1.0_3.1.4 R/4.0.3
  module load gdal geos proj

  # Allele frequency concatenation script
  Rscript ${wd}/3.slim_vcf_output_pairwisefst_pca.R \
  ${tmpdir} \
  ${seed}_.acount \
  ${wd}/data/freq.output.new

  # Move output files
  mv ${tmpdir}/slim_constant-pop* ${output}/
  mv ${wd}/data/freq.output.new* ${output}/

  # Remove VCFs and temporary data
  echo "Removing temporary files"
  rm -fr ${tmpdir}

done

# Remove total tmp directory
rm -fr ${wd}/tmp/${SLURM_JOB_ID}

# Finish
echo "Finish"
