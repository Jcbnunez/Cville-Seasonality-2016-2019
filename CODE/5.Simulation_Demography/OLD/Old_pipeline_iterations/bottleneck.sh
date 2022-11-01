#!/usr/bin/env bash
#
#SBATCH -J slim-constant # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 7G
#SBATCH -o /scratch/csm6hg/slim/slim-constant.out # Standard output
#SBATCH -e /scratch/csm6hg/slim/slim-constant.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load samtools/1.10
module load htslib/1.9
module load intel/18.0 intelmpi/18.0 R/3.6.3
module load gcc/7.1.0 slim

# Working & temp directory
wd="/scratch/csm6hg/slim"
outdir="/scratch/csm6hg/slim/test"

# Extract constants from parameter file
slurmID=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f1 )
EG=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
K=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Rep=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
nSamp=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Gen=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )
seed=$( cat ${wd}/constant_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f7 )

# RAM disk
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

# Run SLiM
slim -d slurmIDi=$slurmID -d EGi=$EG -d Ki=$K -d Repi=$Rep -d nSampi=$nSamp -d Geni=$Gen -d simIDi=$seed -d SJi=${SLURM_JOB_ID} -d STi=${SLURM_ARRAY_TASK_ID} ${wd}/constant.slim

# Gunzip and tabix within subdirectories
for i in ${tmpdir}/*${seed}_.vcf ; do bgzip ${i} ; done
for i in ${tmpdir}/*${seed}_.vcf.gz ; do tabix -f -p vcf ${i} ; done

# Run through PopGenome
echo "Running through PopGenome"

# Diversity script
Rscript ${wd}/slim_vcf_output.R ${tmpdir} ${seed}_.vcf.gz ${outdir}/constant.output.csv

# Population script
Rscript ${wd}/compile.write.popsize.R ${tmpdir} ${seed}_.txt ${outdir}/constant.pop.csv

# Remove vcfs
echo "Removing files"

rm -fr ${tmpdir}

# Finish
echo "Finish"
