#!/usr/bin/env bash
#
#SBATCH -J r2.ag # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/r2.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/r2.%A_%a.err # Standard error

module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

Rscript \
--vanilla \
1.0.calculate_linkage_w_inversion.r

date
echo "done"