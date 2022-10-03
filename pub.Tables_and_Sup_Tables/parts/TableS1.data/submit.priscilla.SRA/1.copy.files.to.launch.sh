#!/usr/bin/env bash
#
#SBATCH -J move.to.sra.port # A single job name for the array
#SBATCH --ntasks-per-node=2 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 15:00:00 ### 
#SBATCH --mem 8G
#SBATCH -o slurmOut/fsts.ag%A_%a.out # Standard output
#SBATCH -e slurmOut/fsts.ag%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#### move files to SRA
#### 
#### 


SRAport=/project/berglandlab/SRA_submission_port

ToSubmit=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/master.sra.file.OW.txt

###


###
#head $r1
#head $r2

###
###
###
###

for i in  {1..85}
do
echo $i

r1=$( cat $ToSubmit  | sed '1d' | awk -F'\t' '{print $7"\t"$8 }' | grep -v "NA" | sed "${i}q;d" | awk '{print $1 }' )
echo $r1

r2=$( cat $ToSubmit  | sed '1d' | awk -F'\t' '{print $7"\t"$8 }' | grep -v "NA" | sed "${i}q;d" | awk '{print $2 }' )
echo $r2

cp $r1 $SRAport
cp $r2 $SRAport 

done





