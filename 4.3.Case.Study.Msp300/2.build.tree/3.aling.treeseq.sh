#### Build tree for Msp300 -- local tree
#### 

#####
#####

module load vcftools
module load samtools
module load tabix
module load bcftools
module load picard
module load tabix
module load bcftools
module load gcc/9.2.0  openmpi/3.1.6
module load mafft
module load iqtree


cat \
Dmel.haps.f/*.fasta \
Dsim.haps.f/*.fasta \
Dsech.XM_032723407.1.fasta \
Dyak.XM_015197606.3.fasta \
> Msp.300.haps.fasta

mafft Msp.300.haps.fasta > Msp.300.haps.al.fasta

mkdir msp.tree
iqtree -s Msp.300.haps.al.fasta -bb 1000 -o D.yak --prefix msp.tree/msp300
