#### Explore Courntey;s VCF
rm(list = ls())

library(foreach)
library(SeqArray)
library(gdsfmt)
library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)

######
tern.gds <- "/project/berglandlab/courtney/simCline/data_files/pooled.gds"
### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen(tern.gds)

snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))

snp.dt %>%
  filter(chr %in% c("Dmel_2L")) %>%
  filter(pos == 5192177 ) %>%
  head


snp.dt %>%
  filter(chr %in% c("Dsim_Scf_2L")) %>% head

