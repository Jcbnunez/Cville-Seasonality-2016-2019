### explore pca CM vs DGRP
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)


load("./dgrp_cm_dimdesc_object.Rdata")

dimdesc_object$Dim.1$quanti %>%
  data.frame() %>% 
  mutate(p.value.adj = p.adjust(p.value, method = "bonferroni"),
         snp_id = rownames(.)) %>% 
  separate(snp_id, into = c("chr", "pos", "feature"), sep = "_" ) %>% 
  ggplot(aes(
    x=as.numeric(pos),
    y=-log10(p.value.adj)
  )) +
  geom_point() ->
  pca_dim1_pval

ggsave(pca_dim1_pval, file = "pca_dim1_pval.pdf")
  
####
dimdesc_object$Dim.1$quanti %>% 
  data.frame() %>% 
  mutate(p.value.adj = p.adjust(p.value, method = "bonferroni"),
         snp_id = rownames(.)) %>% 
  separate(snp_id, into = c("chr", "pos", "feature"), sep = "_", remove = F ) ->
  PCA_dat_parsed

PCA_dat_parsed %>%
  ggplot(aes(
    x=as.numeric(pos),
    y=correlation^2,
    color = -log10(p.value.adj)
  )) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0.85) +
  scale_color_gradient2(low= "firebrick", mid = "purple", high = "steelblue", midpoint = -log10(1e-50) ) +
  theme_bw() +
  theme(legend.position = "top") +
  annotate("text", x = 8e6, y = 0.90, label = "Inv Marker > 0.85") ->
  pca_dim1_cor

ggsave(pca_dim1_cor, file = "pca_dim1_cor.png", w=5, h= 4)


### What are the markers
PCA_dat_parsed %>%
  filter(p.value.adj < 0.01,
         correlation^2 > 0.85) ->
  inversion_markers

###Bring in old markers
svm_2lt_markers = fread("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/CM_2LT_markers/in2lt_ld_informative_markers_FINAL.txt", header = F)

intersect(PCA_dat_parsed$snp_id, svm_2lt_markers$V1)


