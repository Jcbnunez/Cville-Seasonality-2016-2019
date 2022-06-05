### Plot trajectory of haplotypes

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


setwd("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/")
### load haplotype info
load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata")
#sim_polarity,
#windws_snp_matrix_clean,
#joint_figure_polarized_missingdat_clean,

joint_figure_polarized_missingdat_clean %>%
  group_by(SNP_id) %>%
  slice_head() %>% 
  separate(SNP_id, remove = F, into = c("chr", "pos", "type"), sep = "_") %>%
  .$pos %>% as.numeric %>% sort() -> snps_of_interest

#### STOPPED HERE
load("../7.LD/merged.ld.Rdata")
### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B)) 
##
##ld_df %>%
##  as.data.frame() %>%
##  filter(R2 > 0.6) %>%
##  group_by(pair_id) %>%
##  mutate(anchor = min(BP_A, BP_B)) %>% 
##  mutate(pair = case_when(anchor == BP_A ~ BP_B,
##                          anchor != BP_A ~ BP_A)) %>% 
##  mutate(win = case_when(
##    anchor > 5155762 & pair < 5255762 ~ "win5",
##    anchor > 6255762 & pair < 6355762 ~ "win6",
##    anchor > 9505762 & pair < 9605762 ~ "win9")) %>%
##  filter(!is.na(win)) %>%
##  left_join(dplyr::select(sim_polarity, pair = pos, tidy_annot)) %>% 
##  group_by(anchor, win, tidy_annot) %>% 
##  summarise(N = n()) %>% 
##  filter(tidy_annot %in% c("NS", "UTR") ) %>%
##  group_by( win) %>% 
##  slice_max(N, with_ties = F) ->
##  target_SNPS_winld_annot

haplo_windows <- "/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/final.windows.pos.Rdata"
load(haplo_windows)

final.windows.pos


ld_df %>%
  filter(BP_A != BP_B) %>%
  filter(R2 > 0.60) %>%
  filter(BP_diff < 2e5) %>% 
  group_by(pair_id) %>%
  mutate(anchor = min(BP_A, BP_B)) %>% 
  mutate(pair = case_when(anchor == BP_A ~ BP_B,
                          anchor != BP_A ~ BP_A)) %>% 
  mutate(win = 
  case_when(
    anchor > 4650065 & pair < 4799922 ~ "win_4.7",
    anchor > 5100324 & pair < 5349218 ~ "win_5.2",
    anchor > 6100321 & pair < 6349489 ~ "win_6.2",
    anchor > 9500286 & pair < 9700005 ~ "win_9.6"
  )) %>%
  filter(!is.na(win)) %>% 
  left_join(dplyr::select(sim_polarity, pair = pos, tidy_annot)) ->
  haplo_tag_snps_r2
  
haplo_tag_snps_r2 %>% group_by(win) %>% summarize(N= n())

haplo_tag_snps_r2 %>%
  mutate(win = 
           case_when(
             anchor > 4650065 & anchor < 4799922 ~ "win_4.7",
             anchor > 5100324 & anchor < 5349218 ~ "win_5.2",
             anchor > 6100321 & anchor < 6349489 ~ "win_6.2",
             anchor > 9500286 & anchor < 9700005 ~ "win_9.6"
           )) %>%
  group_by(anchor, win) %>%
  summarize(Mean_R2 = mean(R2),
            Npairs = n()) -> haplotag_summaries

haplotag_summaries %>%
  group_by(win) %>% 
  filter(Mean_R2 > 0.6) %>%
  slice_max(Npairs) -> 
  Anchor_snps_of_haplotags

haplo_tag_snps_r2 %>%
  filter(anchor %in% Anchor_snps_of_haplotags$anchor) %>%
  ggplot(
    aes(
      x=pair,
      y=R2,
      color = as.factor(tidy_annot)
    )
  ) + 
  geom_point() +
  geom_hline(yintercept = 0.60) +
  #scale_color_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 0.5) +
  facet_grid(anchor~win, scales = "free_x") ->
  ld_Segments

ggsave(ld_Segments, file = "ld_Segments.pdf", w = 8, h = 8)


### Select haplo-tag snps
haplo_tag_snps_r2 %>%
  filter(anchor %in% Anchor_snps_of_haplotags$anchor) -> all_haplo_pairs

unique(c(
unique(all_haplo_pairs$BP_A),
unique(all_haplo_pairs$BP_B))) %>%
  .[order(.)] -> hap_snps_tag_select

##include breakpoints
joint_figure_polarized_missingdat_clean %>% 
  filter(win %in% c("left_w", "right_w")) %>%
  .$SNP_id %>%
  unique -> inversion_snps_tag_select

### haplo-tag SNPs list
  rbind(
  data.frame(SNP_id = paste("2L", hap_snps_tag_select, "SNP", sep = "_"), class = "GLM_LD"), 
  data.frame(SNP_id = inversion_snps_tag_select, class = "INV")) ->
    haplo_tags_SNPids

save(haplo_tags_SNPids, file = "haplo_tags_SNPids.Rdata")  

##### PLot the haplomap
##joint_figure_polarized_missingdat_clean %>%
##  filter(SNP_id %in% haplo_tags_SNPids$SNP_id) ->
##  joint_figure_polarized_haplotags
##  
########## start plotting
##
##set.seed(123456)
##select_haps = data.frame( 
##  rbind(filter(samp_names, type == "INV"), 
##        sample_n(filter(samp_names, type == "STD" ), 30, replace = FALSE))
##)
##
##joint_figure_polarized_haplotags %>%
##  filter(samp_id %in% select_haps$samp) %>% 
##  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
##  windws_snp_matrix_haplotags
##
##rownames(windws_snp_matrix_haplotags) = windws_snp_matrix_haplotags$samp_id
##
##analyses_types = list(
##  all=c("left_w",
##        "win_4.6",
##        "win_5.1",
##        "win_6.2",
##        "win_6.8",
##        "win_9.5",
##        "right_w"),
##  win4.6=c("win_4.6"),
##  win5.1=c("win_5.1"),
##  win6.2=c("win_6.2"),
##  win6.8=c("win_6.8"),
##  win9.5=c("win_9.5")
##)
##
##foreach(i=1:length(analyses_types))%do%{
##  
##  analysis = names(analyses_types)[i]
##  target_snps <- filter(sim_polarity, win %in% analyses_types[[i]] )$SNP_id 
##  data_in <- windws_snp_matrix_haplotags[, which(colnames(windws_snp_matrix_haplotags) %in% target_snps)  ]
##  actual_sel_snps <- colnames(data_in)
##  rownames(data_in) <- rownames(windws_snp_matrix_haplotags)
##  
##  # dplyr::select(windws_snp_matrix_clean, !(samp_id), filter(SNP_guide_metadata, window == "win5" )$SNP_id ##)
##  D_all <- dist( data_in )
##  tre_all <- njs(D_all)
##  
##  tree_all_plot <- ggtree(tre_all, options(ignore.negative.edge=TRUE)) + 
##    #geom_tiplab(size=2, align=TRUE, linesize=.5) + 
##    theme_tree2() 
##  
##  tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_haplotags + geom_tippoint(aes(color=pop))
##  
##  is_tip <- tre_all$edge[,2] <= length(tre_all$tip.label)
##  ordered_tips <- tre_all$edge[is_tip, 2]
##  tre_all$tip.label[ordered_tips]  -> tree_order
##  
##  joint_figure_polarized_haplotags %>% 
##    filter( SNP_id %in% colnames(windws_snp_matrix_haplotags)[-1] ) %>% 
##    #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
##    mutate(samp_id_fct = factor(samp_id, levels = tree_order)) %>% 
##    ggplot(
##      aes(
##        x=as.factor(loci_id),
##        y=samp_id_fct,
##        fill = polarity
##      )
##    ) + geom_tile(size = 0.1) +
##    facet_grid(.~win, scales = "free", space = "free"
##               #ncol = 1, shrink = F
##    ) +
##    ggtitle(analysis) +
##    theme_bw() +
##    scale_fill_brewer(palette = "Set1") +
##    theme(axis.title.y=element_blank(),
##          axis.text.y=element_blank(),
##          axis.ticks.y=element_blank(),
##          axis.text.x=element_blank(),
##          axis.ticks.x=element_blank())  ->
##    karyo_plot_joint_all
##  
##  karyo_plot_joint_all %>% 
##    insert_left(tree_all_plot, width = 0.1)  -> hap_tree_plots_all
##  
##  ggsave(hap_tree_plots_all, file = paste(analysis, "hap_tree_plots.haplotags.pdf", sep = "."), h = 6, w = ##10)
##  
##  ### add annotation
##  joint_figure_polarized_haplotags %>%
##    filter( SNP_id %in% colnames(windws_snp_matrix_haplotags)[-1] ) %>% 
##    group_by(SNP_id) %>%
##    slice_head() %>%
##    ggplot(aes(
##      x=as.factor(loci_id),
##      y=1,
##      fill = tidy_annot
##    )) + geom_tile(size = 0.1) +
##    facet_grid(.~win, scales = "free", space = "free"
##    ) +
##    #ggtitle(analysis) +
##    theme_bw() +
##    theme(axis.title.y=element_blank(),
##          axis.text.y=element_blank(),
##          axis.ticks.y=element_blank(),
##          axis.text.x=element_blank(),
##          axis.ticks.x=element_blank(),
##          legend.position = "bottom")  ->
##    annots_plot
##  ggsave(annots_plot, file =  "annots_plot.pdf",  h = 1.5, w = 10)
##
##}
##
##
##
##
##
##
##
########## STOPPED HERE
########## STOPPED HERE
########## STOPPED HERE
######## STOPPED HERE

