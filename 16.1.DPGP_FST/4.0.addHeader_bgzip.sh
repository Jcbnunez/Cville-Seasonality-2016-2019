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


header=/project/berglandlab/jcbnunez/emergency_scratch/dpgp3/vcf.header.txt

module load tabix

files=$(ls ./vcfs_diploidized)

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
cat $header \
vcfs_diploidized/${myarr[k]} \
> vcfs_diploidized_header/fixhead.${myarr[k]}

bgzip vcfs_diploidized_header/fixhead.${myarr[k]}
tabix vcfs_diploidized_header/fixhead.${myarr[k]}.gz

date
echo "done"

