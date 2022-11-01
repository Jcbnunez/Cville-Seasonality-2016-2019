#!/usr/bin/env bash
#
#SBATCH -J recomp # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-20

reads_loc=/project/berglandlab/SRA_submission_port/fix_20/pending

module load gcc/7.1.0 gcc/9.2.0
module load bowtie2/2.2.9

#bowtie2-build /project/berglandlab/Dmel_genomic_resources/References/holo_dmel_6.12.fa holo_dmel_6.bti.bt2

ifile=$(cat ./3.files.to.recom.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $ifile

mkdir folder_${SLURM_ARRAY_TASK_ID}

cd folder_${SLURM_ARRAY_TASK_ID}

bowtie2 --very-fast  -U $reads_loc/$ifile.gz \
-x ../holo_dmel_6.bti.bt2 \
--un ./$ifile.un  \
--al ./$ifile.al \
-S ./$ifile.SAM

#ifile=$(cat files.to.recom.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
#echo $ifile

#decompress
#gzip -d $ifile.gz

#spot check
#head $ifile
#tail $ifile

#recompress
#gzip $ifile

date
echo "done"
