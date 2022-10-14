#!/usr/bin/env bash
#
#SBATCH -J remake.BPGP # A single job name for the array
#SBATCH --ntasks-per-node=5 # 5 cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 
#SBATCH --mem 20G
#SBATCH -o slurmOut/glm.fsts.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/glm.fsts.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-197

#cp -r  /scratch/yey2sn/Supergene_paper/2.DPGP/SEQs/dpgp3_sequences/fasta_outs ./

module load gcc/7.1.0
module load openmpi/3.1.4
module load R/4.1.1
module load gdal
module load proj

module load anaconda
source activate bio_space

module load tabix

#cat /project/berglandlab/Dmel_genomic_resources/References/Dmel5/dm3.2L.liner.fa |
#tr [:lower:] [:upper:]  > /project/berglandlab/Dmel_genomic_resources/References/Dmel5/dm3.2L.liner.UPR.fa
REF=/project/berglandlab/Dmel_genomic_resources/References/Dmel5/dm3.2L.liner.UPR.fa

mkdir alns
mkdir vcfs_o

files=$(ls ./fasta_outs)

myarr=()

for i in ${files[@]}
do
   echo $i
   myarr+=("$i")
done
#echo "${#myarr[@]}"

k=${SLURM_ARRAY_TASK_ID}

###
echo ${myarr[k]}

## Merging files
cat $REF \
fasta_outs/${myarr[k]} \
> alns/aligned.${myarr[k]}.fa

#grep ">" aligned.${myarr[k]}.fa

####
cat alns/aligned.${myarr[k]}.fa | bio format --vcf > vcfs_o/${myarr[k]}.vcf
head -n 20 vcfs_o/${myarr[k]}.vcf

#### implement R script to clean the VCF

Rscript \
--vanilla \
3.1.process.vcfs.R \
vcfs_o/${myarr[k]}.vcf


#bgzip vcfs_o/${myarr[k]}.vcf
#tabix vcfs_o/${myarr[k]}.vcf.gz

date
echo "done"

