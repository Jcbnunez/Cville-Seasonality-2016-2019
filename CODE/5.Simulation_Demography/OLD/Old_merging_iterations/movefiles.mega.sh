#!/usr/bin/env bash
#SBATCH -J mega_merge # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH --cpus-per-task=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 #  hours
#SBATCH --mem 5G
#SBATCH -o /scratch/csm6hg/bottleneck/err/mv.out # Standard output
#SBATCH -e /scratch/csm6hg/bottleneck/err/mv.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

find /project/berglandlab/connor/bottleneck/data3 -name '*pca*' -exec \
mv {} /project/berglandlab/connor/bottleneck/data3.pca \;

find /project/berglandlab/connor/bottleneck/data3 -name '*euc*' -exec \
mv {} /project/berglandlab/connor/bottleneck/data3.euc \;

find /project/berglandlab/connor/bottleneck/data3 -name '*fst*' -exec \
mv {} /project/berglandlab/connor/bottleneck/data3.fst \;

find /project/berglandlab/connor/bottleneck/data3 -name '*afvar*' -exec \
mv {} /project/berglandlab/connor/bottleneck/data3.afvar \;
