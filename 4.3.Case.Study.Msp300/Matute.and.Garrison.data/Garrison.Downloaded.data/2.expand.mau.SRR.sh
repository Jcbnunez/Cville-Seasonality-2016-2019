#! /bin/bash

#SBATCH -J srr.downl 
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 7:00:00 
#SBATCH --mem 40G
#SBATCH -o ./slurmout/srr.downl.%A_%a.out # Standard output
#SBATCH -e ./slurmout/srr.downl.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account jcbnunez
#SBATCH --array=1-30

module load sratoolkit

guide=/scratch/yey2sn/Overwintering_ms/msp300.case/mauritania.srr.fetch.txt
root=/scratch/yey2sn/Overwintering_ms/msp300.case/mauritania.SRAs
mkdir ./fasta.files.maur


SRR=`sed '1d' $guide | awk -F " " '{print $2}' | sed -n ${SLURM_ARRAY_TASK_ID}p`
echo $SRR

fasterq-dump \
--split-files $root/$SRR/$SRR.sra \
--outdir ./fasta.files.maur

#gzip ./fasta.files.fld/$SRR.sra.fastq