### Quantify number of SNPS
### 

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(gtools)
library(poolfstat)
library(vroom)

load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")

snp.dt %<>% filter(chr != "X")

snp.dt %>% summarise(N= n())

snp.dt %>% group_by(chr) %>% summarise(N= n())

### -- individual sequencing

system("vcf=/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz;  module load bcftools; bcftools query -f '%CHROM\t%POS\n' $vcf > ind.snps.txt ")

ind.snps <- vroom("ind.snps.txt", col_names = F)
ind.snps %>%
  filter(X1 %in% c("2L","2R","3L","3R")) %>% 
  dim

