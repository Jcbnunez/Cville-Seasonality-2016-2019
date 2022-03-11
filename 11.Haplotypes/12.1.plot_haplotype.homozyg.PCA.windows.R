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
library(car)
library(DescTools)
library(ape)
library(ggtree)
library(aplot)
library(forcats)
registerDoMC(2)

###
setwd("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/")
# load the inversion markers
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
# load outlier GLMs
glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"
load(glm.file)
glm.out %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(mod == "aveTemp+year_factor",
         chr == "2L",
         rnp.clean < 0.05) -> 
  glm.outliers.2L

#outliers_glm %<>% separate(V1, into = c("chr", "pos", "type"))

#load windows data
load("./AF_dat_summ_id.dat.Rdata")
AF_dat_summ_id$win %>% table
SIM_AF %>% head
###

STD_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt", header = F)
STD_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/STD_DGRP_OnlyNames.txt", header = F)
INV_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt", header = F)
INV_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/INV_DGRP_OnlyNames.txt", header = F)

vcf_file <- "./MC.DGRP.merged.readyForImport.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE)
##my_genind <- vcfR2genind(vcf)
#---> AF_dat_summ_id

left = filter(AF_dat_summ_id, win == "left", SNP_id %in% final_in2Lt_markers)
right = filter(AF_dat_summ_id, win == "right", SNP_id %in% final_in2Lt_markers)
win5 = filter(AF_dat_summ_id, win == "win_5", SNP_id %in% glm.outliers.2L$SNP_id)
win6 = filter(AF_dat_summ_id, win == "win_6", SNP_id %in% glm.outliers.2L$SNP_id)
win9 = filter(AF_dat_summ_id, win == "win_9", SNP_id %in% glm.outliers.2L$SNP_id)

markers_left = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% left$SNP_id)], polyThres = 0.00)
markers_right = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% right$SNP_id)], polyThres = 0.00)
markers_win5 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win5$SNP_id)], polyThres = 0.00)
markers_win6 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win6$SNP_id)], polyThres = 0.00)
markers_win9 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% win9$SNP_id)], polyThres = 0.00)

datasets = list(left_w = markers_left,
                right_w = markers_right,
                win5 = markers_win5,
                win6 = markers_win6,
                win9 = markers_win9)

#################################################
##### PLOT JOINNT FIGURE
joint_figure = foreach(i = 1:length(datasets), .combine = "rbind" )%dopar%{
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
    mutate(win =  names(datasets)[i]) -> objt
  
  return(objt)
  
}

##prepare plot
## polarization step  
left$SNP_id %>% data.frame( SNP_id = . , win = "left_w",  loci_id = 1:length(.)) -> a
right$SNP_id %>% data.frame( SNP_id = . , win = "right_w",  loci_id = 1:length(.)) -> b
win5$SNP_id  %>% data.frame( SNP_id = . , win = "win5",  loci_id = 1:length(.)) -> c
win6$SNP_id  %>% data.frame( SNP_id = . , win = "win6",  loci_id = 1:length(.)) -> d
win9$SNP_id  %>% data.frame( SNP_id = . , win = "win9",  loci_id = 1:length(.)) -> e

SIM_AF %<>% mutate(SNP_id = paste(chr,pos,"SNP", sep = "_"))

rbind(a,b,c,d,e) %>%
  left_join(SIM_AF) -> sim_polarity

joint_figure %>% 
  left_join(sim_polarity) %>%  
  mutate(value_pol = case_when(
    af == 0 ~ 1-as.numeric(value),
    af == 1 ~ as.numeric(value))) %>% 
  mutate(polarity = case_when(
    value_pol == 0 ~ "Ancestral",
    value_pol == 1 ~ "Derived")) ->
  joint_figure_polarized

tot_snps = 680
joint_figure_polarized %>%
  group_by(loci_id,value_pol, win) %>%
  summarize(N= (n()/tot_snps)) %>% 
  filter(is.na(value_pol)) %>% 
  mutate(id_id_win = paste(loci_id, win, sep = "_" )) %>%
  filter(N < 0.2) %>%
  .$id_id_win -> id_id_wins

joint_figure_polarized %>%
mutate(id_id_win = paste(loci_id, win, sep = "_" )) %>%
filter(id_id_win %in% id_id_wins) ->
  joint_figure_polarized_missingdat_clean

joint_figure_polarized_missingdat_clean %>%
.[complete.cases(.$value_pol),] %>%
ggplot(
    aes(
      x=as.factor(loci_id),
      y=samp_id,
      fill = polarity
    )
  ) + geom_tile(size = 0.1) +
  facet_grid(karyo+pop~win, scales = "free", 
             #ncol = 1, shrink = F
             ) +
  ggtitle("joint hapblocks") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  ->
  karyo_plot_joint

ggsave(karyo_plot_joint, file =  "karyo_plot_joint.png", w = 12, h = 8)
ggsave(karyo_plot_joint, file =  "karyo_plot_joint.pdf", w = 12, h = 8)

############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
### ONLY INVERSION
samp_names =
  data.frame(rbind(
    cbind(samp=paste(STD_CM$V1, "0", sep = "_"), type = "STD"),
    cbind(samp=paste(STD_CM$V1, "1", sep = "_"), type = "STD"),
    cbind(samp=paste(STD_DGRP$V1, "0", sep = "_"), type = "STD"),
    cbind(samp=paste(STD_DGRP$V1, "1", sep = "_"), type = "STD"),
    cbind(samp=paste(INV_CM$V1, "0", sep = "_"), type = "INV"),
    cbind(samp=paste(INV_CM$V1, "1", sep = "_"), type = "INV"),
    cbind(samp=paste(INV_DGRP$V1, "0", sep = "_"), type = "INV"),
    cbind(samp=paste(INV_DGRP$V1, "1", sep = "_"), type = "INV")
  ))

set.seed(123456)
select_haps = data.frame( 
  rbind(filter(samp_names, type == "INV"), 
        sample_n(filter(samp_names, type == "STD" ), 20, replace = FALSE))
)

joint_figure_polarized_missingdat_clean %>%
  filter(samp_id %in% select_haps$samp) %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_clean

rownames(windws_snp_matrix_clean) = windws_snp_matrix_clean$samp_id

D_all <- dist(dplyr::select(windws_snp_matrix_clean, !(samp_id)))
tre_all <- nj(D_all)

tree_all_plot <- ggtree(tre_all) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

gheatmap(tree_all_plot, dplyr::select(windws_snp_matrix_clean, !(samp_id)), offset=12, width=5.0, 
         colnames=FALSE, legend_title="genotype") +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3)) -> tree_all_genos
ggsave(tree_all_genos, file = "tree_all_genos.pdf", h = 6, w = 10)

#### only CM
#### only CM
#### only CM
#### only CM
#### only CM
#### only CM
#### only CM

joint_figure_polarized_missingdat_clean %>%
  filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_clean_CM

rownames(windws_snp_matrix_clean_CM) = windws_snp_matrix_clean_CM$samp_id

D_cm <- dist(dplyr::select(windws_snp_matrix_clean_CM, !(samp_id)))
tre_cm <- nj(D_cm)

tree_cm <- ggtree(tre_cm) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

is_tip <- tre_cm$edge[,2] <= length(tre_cm$tip.label)
ordered_tips <- tre_cm$edge[is_tip, 2]
tre_cm$tip.label[ordered_tips]  -> tree_order

joint_figure_polarized_missingdat_clean %>%
  filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
  mutate(samp_id_fct = factor(samp_id, levels = tree_order)) %>%
  ggplot(
    aes(
      x=as.factor(loci_id),
      y=samp_id_fct,
      fill = polarity
    )
  ) + geom_tile(size = 0.1) +
  facet_grid(.~win, scales = "free", space = "free"
             #ncol = 1, shrink = F
  ) +
  ggtitle("joint hapblocks") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  ->
  karyo_plot_joint_CM

karyo_plot_joint_CM %>% insert_left(tree_cm, width = 0.1) -> hap_tree_plots

ggsave(hap_tree_plots, file = "hap_tree_plots.cm.pdf", h = 6, w = 10)

### DO PCA
### DO PCA
### DO PCA
### DO PCA
### DO PCA
### DO PCA
### DO PCA
### DO PCA
### DO PCA

joint_figure_polarized_missingdat_clean %>%
  #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_clean_all

rownames(windws_snp_matrix_clean_all) = windws_snp_matrix_clean_all$samp_id

dplyr::select(windws_snp_matrix_clean_all, !(samp_id)) %>%
  PCA(graph = F) ->
  PCA_obj_tmp

PCA_obj_tmp$ind$coord %>%
  as.data.frame() %>%
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
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.2,
      color = as.factor(karyo),
      shape = pop
    )
  ) +
  geom_point(size = 2) ->
  pca_plot

ggsave(pca_plot, file = "all.samps.pca_plot.pdf")

matrix_with_ids %>%
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.3,
      color = as.factor(karyo),
      shape = pop
    )
  ) +
  geom_point(size = 2) ->
  pca_plot3

ggsave(pca_plot3, file = "all.samps.pca_plot3.pdf")

matrix_with_ids %>%
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.4,
      color = as.factor(karyo),
      shape = pop
    )
  ) +
  geom_point(size = 2) ->
  pca_plot4

ggsave(pca_plot4, file = "all.samps.pca_plot4.pdf")


#Contribution
scaleFUN <- function(x) sprintf("%.2f", x)
PCA_obj_tmp$var$cor %>%
  as.data.frame() %>%
  dplyr::select(Dim.1, Dim.2) %>%
  mutate(pos_id = rownames(.)) %>% 
  mutate(pos_bin = RoundTo(as.numeric(pos_id), 1e6, "floor") ) %>% 
  filter(pos_bin %in% c(5e6, 6e6, 9e6)) %>%
  melt(id = c("pos_id","pos_bin") ) %>% 
  ggplot(
    aes(
      x=as.numeric(pos_id)/1e6,
      y=value^2,
      color = variable
    )
  ) + geom_point() + 
  scale_x_continuous(labels=scaleFUN) +
  facet_grid(variable~as.factor(pos_bin),scales = "free_x") ->
  pca_contrib_plot

ggsave(pca_contrib_plot, file = "pca_contrib_plot.pdf", h = 3, w = 6 )
