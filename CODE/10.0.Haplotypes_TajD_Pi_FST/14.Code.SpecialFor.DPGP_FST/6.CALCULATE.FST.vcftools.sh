#! /bin/bash

#SBATCH -J fst_inv_std 
#SBATCH --ntasks-per-node=10 # 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account jcbnunez


# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load vcftools
module load bcftools
module load tabix


#Final VCF file
input_vcf=./DPGP3.2L.merged.vcf.gz

KAR=DPGP3

#File containing populations
std_samps=std.dpgp3.txt
inv_samps=inv.dpgp3.txt


#window parameters
window=100000
step=50000

window2=500000
step2=100000

window3=10000
step3=5000

window4=5000
step4=1000


#######
vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window \
--fst-window-step $step \
--chr CHR2L \
--out FST.W_$window.S_$step.$KAR


vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window2 \
--fst-window-step $step2 \
--chr CHR2L \
--out FST.W_$window2.S_$step2.$KAR

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window3 \
--fst-window-step $step3 \
--chr CHR2L \
--out FST.W_$window3.S_$step3.$KAR

vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window4 \
--fst-window-step $step4 \
--chr CHR2L \
--out FST.W_$window4.S_$step4.$KAR


###
#vcftools \
#--gzvcf $input_vcf \
#--from-bp 5182177 \
#--to-bp 5202177 \
#--chr CHR2L \
#--recode --recode-INFO-all \
#--out msp300

