#!/usr/bin/env bash
#
#SBATCH -J ld_plink # A single job name for the array
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 3:00:00 ### 3 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/bal7cg/score_error/plink.ld.%A_%a.err # Standard error
#SBATCH -e /scratch/bal7cg/score_output/plink.ld.%A_%a.out # Standard output
#SBATCH -p standard
#SBATCH --mail-type=end
#SBATCH --mail-user=bal7cg@virginia.edu
#SBATCH --account berglandlab_standard

### This script output a seres of square LD matrices from the DGRP
### Prepared by Dr. Joaquin C. B. Nunez, University of Virginia, May 2022
### yey2sn@virginia.edu

##########
##########
# README #
##########
##########

## DEPENDENCIES

#### ---> This script depends on PLINK 1.9
#### ---> This script depends on an external R script called "make_square_matrix.R"; distributed also with this code
#### ---> The R dependency requires libraries: data.table(1.14.2), tidyverse(1.3.1), magrittr(2.0.1)

#### ---> NOTICE THAT:
## This script requires a series of input files properly placed in pre-made folders. 
## Failure to prepare these files/folder will result in the script failing. 

## Needed files
##1: the VCF of the DRGP.. or any input VCF for that matter!
### found here:
vcf=/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.vcf.gz
### WARNING ^^^ above is ACTUAL CODE. DO NOT COMMENT IT OUT! ^^^ WARNING ###

##2: A master file with all the job tasks, herein called "./master_job_task.txt"
## This file contains the name of job batches (less than 9999 jobs due to HPC policy) where
## each job covers at least a couple hundred genome windows for LD calculation. 
## This file does not have a header and has two columns; tab separated.
## The first column is the job number. The second column is a file address. Inside each file
## there should be the addresses of various window files
#1	jobs/job1.txt
#2	jobs/job2.txt
#3	jobs/job3.txt
#4	jobs/job4.txt
#5	jobs/job5.txt

## like so:
master_file=./master_job_task.txt
### WARNING ^^^ above is ACTUAL CODE. DO NOT COMMENT IT OUT! ^^^ WARNING ###

##3: A folder called "jobs" where the job files should be stored... e.g. "jobs/job3.txt"

##4: What is inside each job file?
## each of these files just have one column. No header.
## Inside each job file there should contain the addresses to window files. In these context,
## window files will contain the SNPid of mutations inside a window of interest.
#./windows/win.103.txt
#./windows/win.104.txt
#./windows/win.105.txt
#./windows/win.106.txt
#./windows/win.107.txt

##5: A folder called "Windows" populated with window files.
## These window files should have the exact naming conventions of the VCF. Check each VCF
## individually to discover this. If this is new data (i.e., if this is applied to non DGRP data)
##, make sure to annotate the VCF with SNPids... GATK can do this, for example...
## For the DGRP the SNP ids look like this:
#2L_5105071_SNP
#2L_5105137_SNP
#2L_5105142_SNP
#2L_5105148_SNP
#2L_5105154_SNP

##FAQ:
### How many job batches do I need? How many windows? how many SNPs per windows? ...
### A: The short answer is ... it depends. The number and size of windows is analysis dependent
### The number of jobs and windows per job depends on the HPC policy. For UVA;s hpc, 
## One can only launch job arrays of 9999 jobs... so do the math how many windows can you fit in
## 9999 jobs .... etc..

### How is this job supposed to be ran? 
### A: I made this assuming SLURM and an array job ""--array=1-9999" ... across windows.


######################
######################
# Script begins here #
######################
######################

### Step 1
### Load necessary modules
module load plink

### Make the output folder
mkdir output

### Step 2
### extract the window set
window_set=$(cat $master_file | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk '{print $2}' )
echo $window_set

### Step 3
### find the number of windows in the file 
ith_max=$(wc -l  ./$window_set | cut -d ' ' -f 1)
re='^[0-9]+$'
if ! [[ $ith_max =~ $re ]] ; then
   echo "error: Not a number" >&2; 
   else echo "ith_max is numeric -- proceed" ## --> exit 1
fi
ith_max_adj=$(expr $ith_max + 1)
echo $ith_max_adj

#########################
#########################
### check numeric operator
#########################
#########################

### Step 4: launch the loop across windows .

for k in seq $(echo $ith_max_adj)
do
 echo "currently iterating over window" $k "of" $(echo $ith_max_adj)
 window_target=$(cat $window_set | sed "${k}q;d" )  
 echo "now processing" $window_target "which contains" $( wc -l $window_target | cut -d ' ' -f 1) "snps"
 
 jobid=$(echo $window_target | sed 's|./windows/||' | sed 's|.txt||' )
 echo $jobid
  
  	  #use plink to estimate ld among the snp set
   	  plink --vcf $vcf \
  		--allow-extra-chr \
  		--double-id \
  		--r2 \
  		--ld-snp-list $window_target \
  		--ld-window 99999 \
  		--ld-window-r2 0.0 \
  		--out ./output/$jobid

  		#clean up
        rm ./output/$jobid.ld.nosex
        rm ./output/$jobid.ld.log
        
    #apply a conversion script to make plink object into a square matrix
    mkdir square_matrices_final
    
    Rscript \
    --vanilla \
    1.1.make_square_matrix.R \
    ./output/$jobid.ld \
    $jobid \
    $window_target
	
done

echo "im done!" date
