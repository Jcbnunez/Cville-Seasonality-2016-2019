#! /bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

### index genomes

module load gcc/9.2.0
module load bwa/0.7.17

echo "index sechelia"

bwa index GCF_004382195.2_Dsechelia_ASM438219v2_genomic.fna

echo "done"
