#!/usr/bin/env bash
#
#SBATCH -J slim-bottleneck # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 #  hours
#SBATCH --mem 2G
#SBATCH -o /scratch/csm6hg/slim_bottleneck/err/bottleneck.%A_%a.out # Standard output
#SBATCH -e /scratch/csm6hg/slim_bottleneck/err/bottleneck.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Working directory
wd="/scratch/csm6hg/slim_bottleneck"

# Load modules
module load samtools/1.10 parallel
module load htslib/1.9
module load gcc/9.2.0 slim

# Parameter file
paramFile=${wd}/model_paramList5

# Extract constants from parameter file
slurmID=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
nMax=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
nMin=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${paramFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )

# Extract chromosome from each individual
repeatSim () {
slurmID=${1}
nMax=${2}
nMin=${3}
Rep=${4}
nSamp=${5}
Gen=${6}
tmpdir=${7}
replicate=${8}
SJi=${9}
STi=${10}

    # Progress message
    echo "Replicate:" $replicate

    # Load modules
    module load samtools/1.10 parallel
    module load htslib/1.9
    module load gcc/9.2.0 slim

    # Working directory
    wd="/scratch/csm6hg/slim_bottleneck"

    # Output directory
    output="/project/berglandlab/connor/BACKUP_scratch/slim_bottleneck/data"

    # Random seed
    seed=$(($RANDOM))

    # RAM disk set up
    [ ! -d ${wd}/tmp/ ] && mkdir ${wd}/tmp/
    [ ! -d ${wd}/tmp/${SJi} ] && mkdir ${wd}/tmp/${SJi}
    [ ! -d ${wd}/tmp/${SJi}/${STi} ] && mkdir ${wd}/tmp/${SJi}/${STi}
    [ ! -d ${wd}/tmp/${SJi}/${STi}/${seed} ] && mkdir ${wd}/tmp/${SJi}/${STi}/${seed}

    # Temporary directory for VCFs and metadata
    tmpdir=${wd}/tmp/${SJi}/${STi}/${seed}

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
    -d SJi=$SJi \
    -d STi=$STi \
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
    module load goolf/7.1.0_3.1.4 R/4.0.3
    module load gdal geos proj

    # Allele frequency concatenation script
    Rscript ${wd}/3.slim_vcf_output_pairwisefst_pca.R \
    ${tmpdir} \
    ${seed}_.acount \
    ${wd}/data/freq.output

    # Move output files
    mv ${tmpdir}/slim_constant-pop* ${output}/
    mv ${tmpdir}/freq.output* ${output}/

    # Remove VCFs and temporary data
    echo "Removing temporary files"
    rm -fr ${tmpdir}

}

# Export function
export -f repeatSim

# Extract chromosome from every individual
parallel -j 10 repeatSim ::: \
$slurmID ::: \
$nMax ::: \
$nMin ::: \
$Rep ::: \
$nSamp ::: \
$Gen ::: \
$tmpdir ::: \
$( seq 1 ${Rep}) ::: \
${SLURM_JOB_ID} ::: \
${SLURM_ARRAY_TASK_ID}

# Remove total tmp directory
rm -fr ${wd}/tmp/${SLURM_JOB_ID}

# Finish
echo "Finish"
