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

#### Plot outliers only

# load outlier GLMs
outliers <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

#####
load("./AF_dat_summ_id.Rdata")
#---> AF_dat_summ_id
#win5 = filter(AF_dat_summ_id, win == "win_5", pos %in% outliers$pos )
win6 = filter(AF_dat_summ_id, win == "win_6", pos %in% outliers$pos)
win10 = filter(AF_dat_summ_id, win == "win_10", pos %in% outliers$pos)

out_markers_win6n10 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% c(win6$SNP_id,win10$SNP_id) )], polyThres = 0.00)

data_in = out_markers_win6n10

i=1
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
  ggtitle("GLM+LD Outliers in Win6 and Win10 combined") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  ->
  karyo_plot

ggsave(karyo_plot, file = paste( names(datasets)[i] , "karyo_out_6n10_plot.png", sep ="."))



