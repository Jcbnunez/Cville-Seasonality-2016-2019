#!/usr/bin/env bash
#
#SBATCH -J submit.to.sra # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 ### 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o ./sra.%A_%a.out # Standard output
#SBATCH -e ./sra.%A_%a.err # Standard error


module load aspera-connect

ascp -i /home/yey2sn/Aspera_keys/aspera.openssh \
-QT -l100m -k1 -d /project/berglandlab/SRA_submission_port/fix_20 \
subasp@upload.ncbi.nlm.nih.gov:uploads/joaquin.c.b.nunez_gmail.com_v5c4emgx

