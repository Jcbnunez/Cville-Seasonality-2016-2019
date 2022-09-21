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

ToFix=/project/berglandlab/jcbnunez/submit.to.SRA/files.to.fix.txt
###

base=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data/
###
#head $r1
#head $r2

###
###
###
###

for i in  {1..20}
do
echo $i

fix_fol=$( cat $ToFix  | awk -F'\t' '{print $3}' | grep -v "NA" | sed "${i}q;d" )
echo $fix_fol

fix=$( cat $ToFix  | awk -F'\t' '{print $4}' | grep -v "NA" | sed "${i}q;d" )
echo $fix

#r2=$( cat $ToSubmit  | sed '1d' | awk -F'\t' '{print $7"\t"$8 }' | grep -v "NA" | sed "${i}q;d" | awk '{print $2 }' )
#echo $r2

ls $base/$fix_fol/$fix

cp $base/$fix_fol/$fix $SRAport
#cp $r2 $SRAport 

done





