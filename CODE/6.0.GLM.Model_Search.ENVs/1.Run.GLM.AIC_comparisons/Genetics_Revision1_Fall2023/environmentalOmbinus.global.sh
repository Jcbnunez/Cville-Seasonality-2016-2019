#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=8 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 50G
#SBATCH -o /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account biol4559-aob2x

module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

# sbatch --array=1-5000%100 /scratch/aob2x/Overwintering_18_19/EnvironmentalOmnibus_Global/redoPerm/environmentalOmbinus.global.sh
# 1-5000%100
# sacct --format='JobID%30,JobName,Priority,Submit%12,Start%12,Elapsed,NCPU,CPUTime,ExitCode,State' -j 52367172 | grep -v "COMPLETE"
# cat /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_52367172_4822.out
# head /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 10:36:25 EST 2021
# tail /scratch/aob2x/lasso_dest/slurmOut/cm_lasso_28561442_1.out # Fri Nov 12 11:14:06 EST 2021
# sacct -j 52367172 | grep -E "RUNNING|PENDING"
### SLURM_ARRAY_TASK_ID=1
# sacct --format='JobID%30,JobName,Priority,Submit%12,Start%12,Elapsed,NCPU,CPUTime,ExitCode,State' -j 52680484 | grep -v "COMPLETE" | grep -v "batch" | grep "blsmm" | sed 's/ \{1,\}/,/g' | cut -f2 -d',' | cut -f2 -d'_' | tr '\n' ','
# sbatch --array=1609,1617,1640,1641,1642,3334,3338,3371,3920,4467 /scratch/aob2x/Overwintering_18_19/EnvironmentalOmnibus_Global/redoPerm/environmentalOmbinus.global.sh
# sacct -j 52701093


date

cd /scratch/aob2x/

if [ -f /scratch/aob2x/environmental_ombibus_global_permV2/job${SLURM_ARRAY_TASK_ID}.Rdata ]; then
  ls -lh /scratch/aob2x/environmental_ombibus_global_permV2/job${SLURM_ARRAY_TASK_ID}.Rdata
else
  Rscript Overwintering_18_19/EnvironmentalOmnibus_Global/redoPerm/environmentalOmbinus.global.v2.R ${SLURM_ARRAY_TASK_ID} 5000
fi

date
