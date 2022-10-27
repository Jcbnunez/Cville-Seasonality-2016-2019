### 

library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)

### load samp info
homozyg_samps <- fread("./Homozyg_samps_metadat.txt", header=F)
names(homozyg_samps) = c("samp","inv_st")
het_samps <- fread("./Het_samps_metadat.txt", header=F)
names(het_samps) = c("samp","inv_st")
rbind(homozyg_samps, het_samps) -> metadata_tab

dgrp_inv_info =fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
dgrp_inv_info$`DGRP Line` = gsub("DGRP", "line", dgrp_inv_info$`DGRP Line`)
names(dgrp_inv_info)[1] = "samp_hap"

####
################
################v
#### hap net
vcf_cm <- read.vcfR("./CM_pops.2L.invReg.phased.recode.vcf.gz", verbose = TRUE)
vcf_drgp <- read.vcfR("./DGRP.2L.invReg.phased.recode.vcf.gz", verbose = TRUE)

genlight_cm <- vcfR2genlight(vcf_cm)
genlight_dgrp <- vcfR2genlight(vcf_drgp)

#filter
tab(genlight_dgrp) %>% 
  colnames(.) %>%
  .[grep("SNP", .)]-> retained_snps_dgrp
tab(genlight_cm) %>% 
  colnames(.) %>%
  .[grep("SNP", .)]-> retained_snps_dcm

intersect(retained_snps_dgrp,retained_snps_dcm) -> intersect_snps

tab(genlight_cm) %>%
  as.data.frame() %>%
  .[, which(names(.) %in% intersect_snps)] -> genlight_cm_subset

tab(genlight_dgrp) %>%
  as.data.frame() %>%
  .[, which(names(.) %in% intersect_snps)] -> genlight_dgrp_subset

dim(genlight_cm_subset)
dim(genlight_dgrp_subset)

rbind(genlight_cm_subset, genlight_dgrp_subset) ->
  merged_gen_tab

save(merged_gen_tab, file = "merged_genDGRP_CM_tab.Rdata")

#merged_gen_tab[,sample( dim(merged_gen_tab)[2],100)] %>% 
merged_gen_tab %>%
  PCA(., 
      graph = FALSE,
      scale.unit = F,
      ncp=5 ) -> pca_fig_std

save(pca_fig_std, file = "pca_fig_std.DGRP.CM.Rdata")

pca_fig_std$ind$coord %>% 
  as.data.frame %>% 
  mutate(samp_hap =  rownames(.)) %>% 
  mutate(dataset = case_when(
    samp_hap %in% rownames(genlight_cm_subset) ~ "CM",
    samp_hap %in% rownames(tab(genlight_dgrp)) ~ "DGRP",
  )) %>%
  left_join(., dgrp_inv_info) -> pca_coords

pca_coords$In_2L_t[is.na(pca_coords$In_2L_t)] = "not_dgrp"
  
pca_coords %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = dataset,
    shape = In_2L_t
  )) +
  geom_point(alpha = 0.5) -> pca_gg_figure_dgrp

ggsave(pca_gg_figure_dgrp, file = "pca_gg_figure_dgrp.pdf")

