library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)
library(foreach)
library(doMC)
registerDoMC(2)

setwd("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/")

STD_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt", header = F)
STD_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/STD_DGRP_OnlyNames.txt", header = F)
INV_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt", header = F)
INV_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/INV_DGRP_OnlyNames.txt", header = F)

vcf_file <- "./MC.DGRP.merged.readyForImport.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE)
##my_genind <- vcfR2genind(vcf)

load("./AF_dat_summ_id.Rdata")
#---> AF_dat_summ_id
#
#
#
#
left = filter(AF_dat_summ_id, win == "left")
right = filter(AF_dat_summ_id, win == "right")
win5 = filter(AF_dat_summ_id, win == "win_5")
win6 = filter(AF_dat_summ_id, win == "win_6")
win10 = filter(AF_dat_summ_id, win == "win_10")

markers_left = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% left$SNP_id)], polyThres = 0.00)
markers_right = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% right$SNP_id)], polyThres = 0.00)
markers_win5 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win5$SNP_id)], polyThres = 0.00)
markers_win6 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win6$SNP_id)], polyThres = 0.00)
markers_win10 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win10$SNP_id)], polyThres = 0.00)


datasets = list(left_w = markers_left,
                right_w = markers_right,
                win5 = markers_win5,
                win6 = markers_win6,
                win10 = markers_win10)

foreach(i = 1:length(datasets) )%dopar%{
#for(i in 1:length(datasets)){
data_in = datasets[[i]]

message(names(datasets)[i])

data_in@tab %>%
  as.data.frame %>% 
  as.data.frame -> data_in_tab

data_in_tab[,seq(from=1, to=dim(data_in_tab)[2], by=2)] -> data_in_tab_loc
dim(data_in_tab_loc)

colnames(data_in_tab_loc) %>%
  data.frame(samp_hap = .) %>% 
  mutate(loci_id = 1:dim(.)[1]) ->  metadat_loci

data_in_tab_loc %>%
  mutate(samp_id = rownames(.)) %>%
  melt(id = c("samp_id"), variable.name = "loci_id") ->
  data_in_tab_loc_melt

data_in_tab_loc_melt$loci_id = as.integer(data_in_tab_loc_melt$loci_id )

data_in_tab_loc_melt %<>%
  mutate(hap_name = gsub("_0$|_1$", "", data_in_tab_loc_melt$samp_id))

data_in_tab_loc_melt %<>%
  mutate(karyo_samp = case_when(
    hap_name %in% STD_CM$V1 ~ "STD_CM",
    hap_name %in% STD_DGRP$V1 ~ "STD_DGRP",
    hap_name %in% INV_CM$V1 ~ "INV_CM",
    hap_name %in% INV_DGRP$V1 ~ "INV_DGRP"
  )) %>% separate(karyo_samp, into = c("karyo", "pop"), sep = "_")


left_join(data_in_tab_loc_melt, metadat_loci) %>% 
  ggplot(
    aes(
      x=as.numeric(loci_id),
      y=samp_id,
      fill = as.factor(value)
    )
  ) + geom_tile(size = 0.1) +
  facet_wrap(karyo~pop, scales = "free_y", ncol = 1, shrink = F) +
  ggtitle(names(datasets)[i]) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  ->
  karyo_plot

ggsave(karyo_plot, file = paste( names(datasets)[i] , "karyo_plot.png", sep ="."))

}

#### Plot outliers only GLM
#### Plot outliers only GLM
#### Plot outliers only GLM
#### Plot outliers only GLM
#### Plot outliers only GLM

# load outlier GLMs
outliers_glm <- fread("./VA_ch_0.05_ith0_variant_id.txt", head = F)
outliers_glm %<>% separate(V1, into = c("chr", "pos", "type"))

#####
load("./AF_dat_summ_id.Rdata")
#---> AF_dat_summ_id
win5_glm = filter(AF_dat_summ_id, win == "win_5", pos %in% outliers_glm$pos )
win6_glm = filter(AF_dat_summ_id, win == "win_6", pos %in% outliers_glm$pos)
win10_glm = filter(AF_dat_summ_id, win == "win_10", pos %in% outliers_glm$pos)

markers_win5_glm = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win5_glm$SNP_id)], polyThres = 0.00)
markers_win6_glm = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win6_glm$SNP_id)], polyThres = 0.00)
markers_win10_glm = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win10_glm$SNP_id)], polyThres = 0.00)

# load outlier GLMs
outliers_glm_ld <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt", head = T)
#---> AF_dat_summ_id
win5_glm_ld = filter(AF_dat_summ_id, win == "win_5", pos %in% outliers_glm_ld$pos )
win6_glm_ld = filter(AF_dat_summ_id, win == "win_6", pos %in% outliers_glm_ld$pos)
win10_glm_ld = filter(AF_dat_summ_id, win == "win_10", pos %in% outliers_glm_ld$pos)

markers_win5_glm_ld = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win5_glm_ld$SNP_id)], polyThres = 0.00)
markers_win6_glm_ld = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win6_glm_ld$SNP_id)], polyThres = 0.00)
markers_win10_glm_ld = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win10_glm_ld$SNP_id)], polyThres = 0.00)

datasets_glm = list(
  markers_win5_glm = markers_win5_glm,
  markers_win6_glm = markers_win6_glm,
  markers_win10_glm = markers_win10_glm,
  markers_win5_glm_ld = markers_win5_glm_ld,
  markers_win6_glm_ld = markers_win6_glm_ld,
  markers_win10_glm_ld= markers_win10_glm_ld)

#### plot haplotypes
#### plot haplotypes
#### plot haplotypes
foreach(i = 1:length(datasets_glm) )%dopar%{
##foreach(i = 3)%dopar%{
    
  #for(i in 1:length(datasets)){
  data_in = datasets_glm[[i]]
  
  message(names(datasets_glm)[i])
  
  data_in@tab %>%
    as.data.frame %>% 
    as.data.frame -> data_in_tab
  
  data_in_tab[,seq(from=1, to=dim(data_in_tab)[2], by=2)] -> data_in_tab_loc
  dim(data_in_tab_loc)
  
  colnames(data_in_tab_loc) %>%
    data.frame(samp_hap = .) %>% 
    mutate(loci_id = 1:dim(.)[1]) ->  metadat_loci
  
  data_in_tab_loc %>%
    mutate(samp_id = rownames(.)) %>%
    melt(id = c("samp_id"), variable.name = "loci_id") ->
    data_in_tab_loc_melt
  
  data_in_tab_loc_melt$loci_id = as.integer(data_in_tab_loc_melt$loci_id )
  
  data_in_tab_loc_melt %<>%
    mutate(hap_name = gsub("_0$|_1$", "", data_in_tab_loc_melt$samp_id))
  
  data_in_tab_loc_melt %<>%
    mutate(karyo_samp = case_when(
      hap_name %in% STD_CM$V1 ~ "STD_CM",
      hap_name %in% STD_DGRP$V1 ~ "STD_DGRP",
      hap_name %in% INV_CM$V1 ~ "INV_CM",
      hap_name %in% INV_DGRP$V1 ~ "INV_DGRP"
    )) %>% separate(karyo_samp, into = c("karyo", "pop"), sep = "_")
  
  
  left_join(data_in_tab_loc_melt, metadat_loci) %>% 
    ggplot(
      aes(
        x=as.numeric(loci_id),
        y=samp_id,
        fill = as.factor(value)
      )
    ) + geom_tile(size = 0.1) +
    facet_wrap(karyo~pop, scales = "free_y", ncol = 1, shrink = F) +
    ggtitle(names(datasets_glm)[i]) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())  ->
    karyo_plot
  
  ggsave(karyo_plot, file = paste( names(datasets_glm)[i] , "karyo_plot.png", sep ="."))
  
}

#####

### DO PCA

foreach(i = 1:length(datasets_glm) )%dopar%{
  data_in = datasets_glm[[i]]
  
  message(names(datasets_glm)[i])
  samp_processed = names(datasets_glm)[i]
  
  data_in@tab %>%
    as.data.frame %>% 
    as.data.frame -> data_in_tab
  
  data_in_tab[,seq(from=1, to=dim(data_in_tab)[2], by=2)] -> data_in_tab_loc
  dim(data_in_tab_loc)
  
  data_in_tab_loc %>%
    as.data.frame %>%
    mutate(samp_id = rownames(.)) %>%
    mutate(hap_name = gsub("_0$|_1$", "", .$samp_id)) %>%
    mutate(karyo_samp = case_when(
      hap_name %in% STD_CM$V1 ~ "STD_CM",
      hap_name %in% STD_DGRP$V1 ~ "STD_DGRP",
      hap_name %in% INV_CM$V1 ~ "INV_CM",
      hap_name %in% INV_DGRP$V1 ~ "INV_DGRP"
    )) %>% separate(karyo_samp, into = c("karyo", "pop"), sep = "_") ->
    matrix_with_ids
  
  
  matrix_with_ids %>%
    filter(karyo == "INV") %>%
    dplyr::select(-c("samp_id", "hap_name", "karyo","pop") ) %>%
    PCA(graph = F) ->
    PCA_obj_tmp
  
  ### clustering analysis
  set.seed(123)
  km.res <- kmeans( data.frame(PCA_obj_tmp$ind$coord) , 3, nstart = 25)
  res <- hcut(data.frame(PCA_obj_tmp$ind$coord), k = 3, stand = TRUE)
  # Visualize
  tree_cut <-
  fviz_dend(res, rect = TRUE, cex = 0.5,
            k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))
  
  ggsave(tree_cut, file = paste(samp_processed, "tree_cut.pdf", sep ="."))
  
  # 3. Visuaze
  data.frame(km.res$cluster) %>%
    mutate(samp_id = rownames(.)) -> km.samps
  
  PCA_obj_tmp$ind$coord %>%
    as.data.frame %>%
    mutate(samp_id = rownames(.)) %>%
    left_join(matrix_with_ids) %>%
    left_join(km.samps) ->
    obj_metdat
  
  save(obj_metdat, file = paste(samp_processed, "cluster_pca.Rdata", sep ="."))
  
  obj_metdat %>%
    ggplot(
      aes(
        x=Dim.1,
        y=Dim.2,
        color = as.factor(km.res.cluster),
        shape = pop
      )
    ) +
    geom_point(size = 2) ->
    pca_plot
  
  ggsave(pca_plot, file = paste(samp_processed, "pca_plot.pdf", sep ="."))

}


###
###
### ---> invetigate clusters
load("markers_win6_glm.cluster_pca.Rdata")
win6_out = obj_metdat

win6_out %>%
  filter(km.res.cluster == 1) %>%
  .$hap_name %>% unique -> Cluster1_win6
write.table(Cluster1_win6, file = "Cluster1_win6.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

win6_out %>%
  filter(km.res.cluster == 2) %>%
  .$hap_name %>% unique -> Cluster2_win6
write.table(Cluster2_win6, file = "Cluster2_win6.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

win6_out %>%
  filter(km.res.cluster == 3) %>%
  .$hap_name %>% unique -> Cluster3_win6
write.table(Cluster3_win6, file = "Cluster3_win6.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)



