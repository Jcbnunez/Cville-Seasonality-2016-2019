library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)

vcf_file <- "./CM_pops.2L.invReg.phased.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

## Make DNAbin

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE)
my_dnabin1

glm_ld_markers_gen = DNAbin2genind(x = my_dnabin1, polyThres = 0.05)

glm_ld_markers_gen@tab %>%
  as.data.frame %>% 
  as.data.frame -> snp_dat_table_all

snp_dat_table_all[,seq(from=1, to=dim(snp_dat_table_all)[2], by=2)] -> snp_dat_table_loc
dim(snp_dat_table_loc)
names(snp_dat_table_loc) %>% head
names(snp_dat_table_loc) %>% tail
####
####
## bring in the data names
std_haps <- fread("./Std_samps_OnlyNames.txt", header = F)
inv_haps <- fread("./Inv_samps_OnlyNames.txt", header = F)


snp_dat_table_loc %>%
mutate(samp_hap =  rownames(.)) %>% 
  mutate(samp = gsub("_0$|_1$", "", samp_hap)) -> 
  snp_dat_table_loc_named
###  
snp_dat_table_loc_named %>%
  mutate(karyo =case_when(samp %in% std_haps$V1 ~ "std_hom",
                          samp %in% inv_haps$V1 ~ "inv_hom") ) %>%
  filter(karyo %in% c("std_hom", "inv_hom") ) -> snp_dat_table_loc_homkar

#save(snp_dat_table_loc_homkar, file = "snp_dat_table_loc_homkar.Rdata")
load("./snp_dat_table_loc_homkar.Rdata")

length(names(snp_dat_table_loc_homkar)) -> L
L.adj = L-3
names(snp_dat_table_loc_homkar) %>% .[1:L.adj] %>% 
  data.frame(variable = .) %>% 
  separate(variable, into = c("pos", "base"), sep = "\\.", remove = F ) %>%
  mutate(pos_id = 1:dim(.)[1]) -> 
  metadat_pos

rownames(snp_dat_table_loc_homkar) %>%
  data.frame(samp_hap = .) %>% 
  mutate(samp_id = 1:dim(.)[1]) ->  metadat_samps

  
snp_dat_table_loc_homkar %>%
  melt(id = c("samp_hap", "samp", "karyo")) ->
  snp_dat_table_loc_homkar_melt

left_join(snp_dat_table_loc_homkar_melt, metadat_pos) -> 
  snp_dat_table_loc_homkar_melt_pos
left_join(snp_dat_table_loc_homkar_melt_pos, metadat_samps) -> 
  snp_dat_table_loc_homkar_melt_pos_samps

save(snp_dat_table_loc_homkar_melt_pos_samps, file = "snp_dat_table_loc_homkar_melt_pos_samps.Rdata")


###
###
###