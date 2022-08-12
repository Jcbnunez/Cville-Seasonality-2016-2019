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

echo "index yakuba"

bwa index GCF_004382145.1_Dmauritiana_ASM438214v1_genomic.fna

echo "done"
