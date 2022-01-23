#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez
#SBATCH --array=1-7

###########################################################################
###########################################################################
# Load R RIVANNA modules
###########################################################################
###########################################################################

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1

###########################################################################
###########################################################################
# Load some necessities for the array to work
###########################################################################
###########################################################################

pop_names=GML_guide_pops.txt

#extract pop name from guide file
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $pop_names`

###########################################################################
###########################################################################
# Run R
###########################################################################
###########################################################################

Rscript Run_YT_model_pvalue_extractor.R $i