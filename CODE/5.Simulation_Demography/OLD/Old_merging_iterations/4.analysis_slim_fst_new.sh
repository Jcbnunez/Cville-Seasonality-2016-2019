#!/usr/bin/env bash
#SBATCH -J vcf2genlight
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 100G
#SBATCH -o /project/berglandlab/connor/err/bindoutput.out # Standard output
#SBATCH -e /project/berglandlab/connor/err/bindoutput.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# Combine all output from Drosophila simulations

# Modules to load
module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj

# Run script
Rscript /project/berglandlab/connor/scripts/4.analysis_slim_fst_new.R
