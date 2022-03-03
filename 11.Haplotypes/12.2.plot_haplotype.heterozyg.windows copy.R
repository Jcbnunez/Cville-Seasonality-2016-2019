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

vcf_file <- "./MC.DGRP.het.merged.readyForImport.recode.vcf.gz"
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
  mutate(karyo_samp = NA)

data_in_tab_loc_melt$karyo_samp[grep("CM", data_in_tab_loc_melt$hap_name)] = "HET_CM"
data_in_tab_loc_melt$karyo_samp[grep("line", data_in_tab_loc_melt$hap_name)] = "HET_DGRP"

#data_in_tab_loc_melt %<>%
#  mutate(karyo_samp = case_when(
#    hap_name %in% STD_CM$V1 ~ "HET_CM",
#    hap_name %in% INV_DGRP$V1 ~ "HET_DGRP"
#  )) %>% separate(karyo_samp, into = c("karyo", "pop"), sep = "_")

data_in_tab_loc_melt %<>%
  separate(karyo_samp, into = c("karyo", "pop"), sep = "_")

left_join(data_in_tab_loc_melt, metadat_loci) %>% 
  ggplot(
    aes(
      x=as.numeric(loci_id),
      y=samp_id,
      fill = as.factor(value)
    )
  ) + geom_tile(size = 0.1) +
  facet_wrap(karyo~pop, scales = "free_y", ncol = 1, shrink = F) +
  ggtitle(paste(names(datasets)[i], "Het Samps", sep = " " ) ) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  ->
  karyo_plot

ggsave(karyo_plot, file = paste( names(datasets)[i] , "karyo_plot.HET.png", sep ="."))

}




