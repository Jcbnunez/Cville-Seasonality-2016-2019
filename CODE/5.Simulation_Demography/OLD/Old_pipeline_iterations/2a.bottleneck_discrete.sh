#!/usr/bin/env bash
#
#SBATCH -J slim-bottleneck # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:30 # 30 mins
#SBATCH --mem 50G
#SBATCH -o /scratch/csm6hg/slim_bottleneck/err/bottleneck.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/slim_bottleneck/err/bottleneck.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Load modules
module load samtools/1.10
module load htslib/1.9
module load gcc/9.2.0 slim

# Working directory
wd="/scratch/csm6hg/slim_bottleneck"

# Parameter file
paramFile=${wd}/model_paramList4

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
nMax=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
nMin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )
seed=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f7 )
replicate=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f8 )

# RAM disk set up
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

# Temporary directory for VCFs and metadata
tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

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
-d SJi=${SLURM_JOB_ID} \
-d STi=${SLURM_ARRAY_TASK_ID} \
-d replicatei=${replicate} \
${wd}/1a.bottleneck_discrete.slim

# Load plink
module load plink/2.00a20210104

# Creates list for all VCFs
ls -d ${tmpdir}/*vcf > \
${tmpdir}/vcf.list

# Extract allele frequency information
while read i; do {

  # Progress message
  echo "Processing:" ${i}

  # Fix output name
  out=$( echo ${i} | sed 's/\.vcf/\n/g' )

  # Run plink Allele counts
  plink2 --vcf ${i} \
  --freq counts cols=+pos \
  --out ${out}

}
done < ${tmpdir}/vcf.list

# Run through Rscript
echo "Running through summary scripts"

# Load R module
module purge
module load goolf/7.1.0_3.1.4 R/4.0.3
module load gdal geos proj

# Allele frequency concatenation script
Rscript ${wd}/3.slim_vcf_output_fst_pca.R \
${tmpdir} \
${seed}_.acount \
${wd}/data/freq.output

# Move population file
mv ${tmpdir}/*pop* ${wd}/data/.

# Remove VCFs and temporary data
echo "Removing temporary files"
rm -fr ${tmpdir}

# Finish
echo "Finish"
