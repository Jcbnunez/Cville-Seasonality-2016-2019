#####
#####

module load vcftools
module load samtools
module load tabix
module load bcftools
module load picard
module load tabix
module load bcftools


# ----> get data
wget https://zenodo.org/record/154261/files/simulans_multisamp_all_chr.vcf.zip
cp /project/berglandlab/courtney/simCline/refgenomes/simulans/dsim-mod.fasta ./

### bgzipping D sim genome
#bgzip simulans_multisamp_all_chr.vcf
#tabix simulans_multisamp_all_chr.vcf.gz

## slice CM
vcftools \
--gzvcf simulans_multisamp_all_chr.vcf.gz \
--recode \
--recode-INFO-all \
--maf 0.01 \
--out Dsim_all_2L.maf0.01 \
--chr 2L 

bgzip Dsim_all_2L.maf0.01.recode.vcf
tabix Dsim_all_2L.maf0.01.recode.vcf.gz
