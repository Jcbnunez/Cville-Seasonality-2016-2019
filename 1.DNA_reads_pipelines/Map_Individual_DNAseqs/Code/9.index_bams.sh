#!/usr/bin/env bash
#
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=45G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

module load samtools

###########################################################################
#Parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################


guide_file=/scratch/yey2sn/phasing_droso_data/Whatshap/samples_to_index.txt

bamfile=`awk -F "\t" '{print $2}' $guide_file | sed -n ${SLURM_ARRAY_TASK_ID}p`

samtools \
index \
-@ $CPU \
$bamfile


