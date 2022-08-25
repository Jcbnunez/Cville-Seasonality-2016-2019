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

out=/scratch/yey2sn/Overwintering_ms/msp300.case/mauritania.SRAs
guide=/scratch/yey2sn/Overwintering_ms/msp300.case/mauritania.srr.fetch.txt

mkdir $out

SRR=$( cat $guide | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $SRR

prefetch $SRR --output-directory $out

