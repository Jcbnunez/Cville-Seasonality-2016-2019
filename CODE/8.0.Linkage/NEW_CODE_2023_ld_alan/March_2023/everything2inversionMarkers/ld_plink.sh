#!/bin/bash
#
#SBATCH -J plink # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/ld/logs/ld.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/ld/logs/ld.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab




### sbatch /home/aob2x/Overwintering_18_19/ld_alan/ld_plink.sh
### sacct -j 49751805
### cat /scratch/aob2x/Overwintering_18_19/ld_alan/logs/ld.49751619_17*.err
### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err






# ijob -A berglandlab_standard -c20 -p standard --mem=30G
module load plink/1.90b6.16

# cp CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz /scratch/aob2x/ld

cd /scratch/aob2x/ld


plink \
--r2 inter-chr with-freqs yes-really \
--double-id --allow-extra-chr \
--ld-window-r2 0.01 \
--vcf /scratch/aob2x/ld/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz \
--ld-snp-list /project/berglandlab/jcbnunez/Shared_w_Alan/in2lt_ld_47snps_informative_markers.txt \
