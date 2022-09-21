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
library(SeqArray)
library(doMC)
library(car)
library(DescTools)
library(ape)
library(ggtree)
library(aplot)
library(forcats)
registerDoMC(2)

###
setwd("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/")

load("./joint_figure_polarized.Rdata")

### find lines of interest
load("./snp.dosage.Rdata")
listdata %>% head

### load inversion info
inv.stat = fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")

inv.stat %<>%
  mutate(ral_id = gsub("DGRP", "line", inv.stat$`DGRP Line` ) )

######

selected_dgrp_lines = 
foreach(i = c(1,3,5,6), .combine = "rbind")%do%{
  print(i)
  
  win.of.int = names(listdata)[i]
  print(win.of.int)
  
  listdata[[which(names(listdata) == win.of.int)]] %>%
    dplyr::select(!c(ral_id,vars,peak) ) %>% 
    rowMeans(na.rm = T) -> row.means.dosage
  
  listdata[[which(names(listdata) == win.of.int)]] %>%
    dplyr::select(c(ral_id,vars,peak) ) %>%
    mutate(mean.dosage = row.means.dosage) %>%
    left_join(inv.stat) -> win.dosage
  
  win.dosage %<>%
    filter(ral_id != "line_630")
  
  win.dosage %>% 
    filter(In_2L_t %in% c("ST")) %>%
    arrange(mean.dosage) %>%
    slice_max(n = 3, mean.dosage) ->
    ST_select
  
  win.dosage %>% 
    filter(In_2L_t %in% c("INV")) %>%
    arrange(mean.dosage) %>%
    slice_min(n = 3, mean.dosage) ->
    INV_select
  
  rbind( 
    mutate(ST_select, line_info = paste("ST.select", win.of.int, sep = ".")),
    mutate(INV_select, line_info = paste("INV.select", win.of.int, sep = "."))
           )

  #win.dosage %>% 
  #  ggplot(aes(
  #    x= fct_reorder(ral_id, mean.dosage),
  #    y=mean.dosage,
  #    color=In_2L_t
  #  )) +
  #  ggtitle(win.of.int) +
  #  geom_point() +
  #  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5),
  #  ) ->
  #  dosage.pts
  #
  #ggsave(dosage.pts, file = paste(win.of.int, "dosage.pts.pdf", sep = ""), 
  #       w = 12, h = 5)
}


selected_dgrp_lines %>%
  filter(In_2L_t == "ST") %>%
  group_by(peak) %>%
  slice_head() -> ids.std

selected_dgrp_lines %>%
  filter(In_2L_t == "INV") ->
  ids.inv

rbind(ids.std, ids.inv)$ral_id %>% unique() -> ids
  
### filter by num of sites
select_pops = c("CM","PA","ME",
                "DGRP")


#Cameroon       Carib          CM        DGRP      France      Guinea 
#20232        8816      388596      499818       22780       12642 
#ME Netherlands          PA      Zambia 
#190162       37608      265974       72936 
##
##joint_figure_polarized %>%
##  filter(pop %in% select_pops) %>% 
##  group_by(SNP_id) %>% 
##  summarize(ALL= n()) -> allSNPS
##
##joint_figure_polarized %>%
##  filter(is.na(value_pol), pop %in% select_pops) %>% 
##  group_by(SNP_id) %>% 
##  summarize(Miss= n()) -> missSNPS
##
##full_join(allSNPS, missSNPS) -> miss_analysis
##
##miss_analysis$Miss[is.na(miss_analysis$Miss)] = 0 
##
##miss_analysis %<>%
##  mutate(perc_miss = as.numeric(Miss)/as.numeric(ALL)) 
##
##miss_analysis$ALL %>% quantile(0.7) -> tresh_all
##miss_analysis$perc_miss %>% quantile(0.7) -> thres_miss
##
##miss_analysis %>%
##  filter(ALL >= tresh_all,
##         perc_miss <= thres_miss) %>% 
##  .$SNP_id -> id_id_wins
##
##### filter by inds
##### 
##
##joint_figure_polarized %>%
##  filter(pop %in% select_pops) %>% 
##  group_by(hap_name) %>% 
##  summarize(ALL= n()) -> inds_snps_samps_density
##
##inds_snps_samps_density %>%
##  filter(ALL >= 2500) %>% 
##  .$hap_name -> keep_samps
##
##
#######
####### Make filterd object
####### 


 ids = c("line_853", #STD1
           "line_634", #STD2
           "line_189", #STD3
           "line_57", #STD4
           
           "line_348", #inv1 - derived 5 and derived 9
           "line_837", #inv2 - ancestral 5 and derived 9
           "line_748", #inv3 - ancestral 5 and derived 9
           "line_161", #inv4 - derived 5 and ancestral 9
           "line_32", #inv5 - derived 5 and ancestral 9
           "line_386" #inv6 - derived 5 and ancestral 9
         ) 

  
filter(joint_figure_polarized, pop %in% c("CM" #,
                                          #"PA",
                                          #"ME"
                                          ) & karyo == "INV" )$hap_name %>% unique() -> inv.nats

joint_figure_polarized %>%
  filter(pop %in% select_pops) %>%
  filter(hap_name %in% c( ids,inv.nats )  ) ->
  joint_figure_polarized_missingdat_clean

joint_figure_polarized_missingdat_clean %>% head

#######

############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
#
#samp_names =
#  data.frame(rbind(
#    cbind(samp=paste(filter(samps_info_pop, karyo == "STD", 
#                            hap_name %in% keep_samps,
#                            pop %in% select_pops)$hap_name, 
#                     "0", sep = "_"), type = "STD"),
#    cbind(samp=paste(filter(samps_info_pop, karyo == "INV", 
#                            hap_name %in% keep_samps,
#                            pop %in% select_pops)$hap_name, 
#                     "0", sep = "_"), type = "INV"),
#    cbind(samp=paste(filter(samps_info_pop, karyo == "STD", 
#                            hap_name %in% keep_samps,
#                            pop %in% select_pops)$hap_name, 
#                     "1", sep = "_"), type = "STD"),
#    cbind(samp=paste(filter(samps_info_pop, karyo == "INV", 
#                            hap_name %in% keep_samps,
#                            pop %in% select_pops)$hap_name, 
#                     "1", sep = "_"), type = "INV")
#  ))
#
#
#set.seed(123456)
#
#select_haps = data.frame( 
#  rbind(filter(samp_names, type == "INV"), 
#        sample_n(filter(samp_names, type == "STD" ), 30, replace = FALSE))
#)

joint_figure_polarized_missingdat_clean %>%
  .[grep("_0", joint_figure_polarized_missingdat_clean$samp_id ),] %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_clean

### Parts of Figure 4
#save(
#  sim_polarity,
#  samp_names,
#  windws_snp_matrix_clean,
#  joint_figure_polarized_missingdat_clean,
#  file = "/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata"
#)

#load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata")

## add rownames
rownames(windws_snp_matrix_clean) = windws_snp_matrix_clean$samp_id

analyses_types = list(
  all.adam=c("left",
        #"win_4.7",
        "win_5.2",
        #"win_6.2",
        "win_9.6",
        "right"),
  #win4.7=c("win_4.7"),
  win5.2.adam=c("win_5.2"),
  #win6.2=c("win_6.2"),
  win9.6.adam=c("win_9.6")
)

foreach(i=1:length(analyses_types))%do%{
  
  message(names(analyses_types)[i])
  
  analysis = names(analyses_types)[i]
  target_snps <- filter(sim_polarity, win %in% as.character(analyses_types[[i]]) )$SNP_id 
  
  data_in <- windws_snp_matrix_clean[, which(colnames(windws_snp_matrix_clean) %in% target_snps)  ] 
  
  actual_sel_snps <- colnames(data_in)
  
  # dplyr::select(windws_snp_matrix_clean, !(samp_id), filter(SNP_guide_metadata, window == "win5" )$SNP_id )
  D_all <- dist( data_in )
  tre_all <- njs(D_all)
  
  tree_all_plot <- ggtree(tre_all, ignore.negative.edge=TRUE) + 
    #geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    theme_tree2()
  
  tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_missingdat_clean + geom_tippoint(aes(color=pop)) + geom_tiplab(as_ylab=TRUE, color='firebrick', size = 7)
  
  
  is_tip <- tre_all$edge[,2] <= length(tre_all$tip.label)
  ordered_tips <- tre_all$edge[is_tip, 2]
  tre_all$tip.label[ordered_tips]  -> tree_order
  
  joint_figure_polarized_missingdat_clean %>%
    #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
    mutate(samp_id_fct = factor(samp_id, levels = tree_order)) %>% 
    .[complete.cases(.$samp_id_fct),] %>%
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
    ggtitle(analysis) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())  ->
    karyo_plot_joint_all
  
  karyo_plot_joint_all %>% 
    insert_left(tree_all_plot, width = 0.1)  -> hap_tree_plots_all
  
  ggsave(hap_tree_plots_all, file = paste(analysis, "hap_tree_plots.all.pdf", sep = "."), h = 6, w = 10)
  
}
