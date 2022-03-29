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
load("/scratch/yey2sn/Overwintering_ms/Figure4/fig4_pt1.Rdata")
#sim_polarity,
#windws_snp_matrix_clean,
#joint_figure_polarized_missingdat_clean,

joint_figure_polarized_missingdat_clean %>%
  group_by(SNP_id) %>%
  slice_head() %>%
  .$pos -> snps_of_interest

load("../7.LD/merged.ld.Rdata")
### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B)) 

ld_df %>%
  as.data.frame() %>%
  filter(R2 > 0.6) %>%
  group_by(pair_id) %>%
  mutate(anchor = min(BP_A, BP_B)) %>% 
  mutate(pair = case_when(anchor == BP_A ~ BP_B,
                          anchor != BP_A ~ BP_A)) %>% 
  mutate(win = case_when(
    anchor > 5155762 & pair < 5255762 ~ "win5",
    anchor > 6255762 & pair < 6355762 ~ "win6",
    anchor > 9505762 & pair < 9605762 ~ "win9")) %>%
  filter(!is.na(win)) %>%
  left_join(dplyr::select(sim_polarity, pair = pos, tidy_annot)) %>% 
  group_by(anchor, win, tidy_annot) %>% 
  summarise(N = n()) %>% 
  filter(tidy_annot %in% c("NS", "UTR") ) %>%
  group_by( win) %>% 
  slice_max(N, with_ties = F) ->
  target_SNPS_winld_annot

ld_df %>%
  filter(R2 > 0.60) %>%
  filter(BP_diff < 1e5) %>% 
  group_by(pair_id) %>%
  mutate(anchor = min(BP_A, BP_B)) %>% 
  mutate(pair = case_when(anchor == BP_A ~ BP_B,
                          anchor != BP_A ~ BP_A)) %>% 
  mutate(win = case_when(
    anchor > 5155762 & pair < 5255762 ~ "win5",
    anchor > 6255762 & pair < 6355762 ~ "win6",
    anchor > 9505762 & pair < 9605762 ~ "win9")) %>%
  filter(!is.na(win)) %>% 
  left_join(dplyr::select(sim_polarity, pair = pos, tidy_annot)) %>% 
  filter(anchor %in% target_SNPS_winld_annot$anchor ) ->
  haplo_tag_snps_r2
  
haplo_tag_snps_r2 %>% group_by(win) %>% summarize(N= n())

haplo_tag_snps_r2 %>%
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
  facet_wrap(anchor~win, scales = "free", ncol = 1) ->
  ld_Segments

ggsave(ld_Segments, file = "ld_Segments.pdf", w = 8, h = 5)


### Select haplo-tag snps
unique(c(
unique(haplo_tag_snps_r2$BP_A),
unique(haplo_tag_snps_r2$BP_B))) %>%
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

### PLot the haplomap
joint_figure_polarized_missingdat_clean %>%
  filter(SNP_id %in% haplo_tags_SNPids$SNP_id) ->
  joint_figure_polarized_haplotags
  
######## start plotting

set.seed(123456)
select_haps = data.frame( 
  rbind(filter(samp_names, type == "INV"), 
        sample_n(filter(samp_names, type == "STD" ), 20, replace = FALSE))
)

joint_figure_polarized_haplotags %>%
  filter(samp_id %in% select_haps$samp) %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_haplotags

rownames(windws_snp_matrix_haplotags) = windws_snp_matrix_haplotags$samp_id

analyses_types = list(
  all=c("win5","win6","win9","left_w", "right_w"),
  win5=c("win5"),
  win6=c("win6"),
  win9=c("win9"))

foreach(i=1:4)%do%{
  
  analysis = names(analyses_types)[i]
  target_snps <- filter(sim_polarity, win %in% analyses_types[[i]] )$SNP_id 
  data_in <- windws_snp_matrix_haplotags[, which(colnames(windws_snp_matrix_haplotags) %in% target_snps)  ]
  actual_sel_snps <- colnames(data_in)
  rownames(data_in) <- rownames(windws_snp_matrix_haplotags)
  
  # dplyr::select(windws_snp_matrix_clean, !(samp_id), filter(SNP_guide_metadata, window == "win5" )$SNP_id )
  D_all <- dist( data_in )
  tre_all <- njs(D_all)
  
  tree_all_plot <- ggtree(tre_all, options(ignore.negative.edge=TRUE)) + 
    #geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    theme_tree2() 
  
  tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_haplotags + geom_tippoint(aes(color=pop))
  
  is_tip <- tre_all$edge[,2] <= length(tre_all$tip.label)
  ordered_tips <- tre_all$edge[is_tip, 2]
  tre_all$tip.label[ordered_tips]  -> tree_order
  
  joint_figure_polarized_haplotags %>% 
    filter( SNP_id %in% colnames(windws_snp_matrix_haplotags)[-1] ) %>% 
    #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
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
  
  ggsave(hap_tree_plots_all, file = paste(analysis, "hap_tree_plots.haplotags.pdf", sep = "."), h = 6, w = 10)
  
  ### add annotation
  joint_figure_polarized_haplotags %>%
    filter( SNP_id %in% colnames(windws_snp_matrix_haplotags)[-1] ) %>% 
    group_by(SNP_id) %>%
    slice_head() %>%
    ggplot(aes(
      x=as.factor(loci_id),
      y=1,
      fill = tidy_annot
    )) + geom_tile(size = 0.1) +
    facet_grid(.~win, scales = "free", space = "free"
    ) +
    #ggtitle(analysis) +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom")  ->
    annots_plot
  ggsave(annots_plot, file =  "annots_plot.pdf",  h = 1.5, w = 10)

}







######## STOPPED HERE
######## STOPPED HERE
######## STOPPED HERE
######## STOPPED HERE

rownames(windws_snp_matrix_clean) -> samps_vec
colnames(windws_snp_matrix_clean) -> snps_vec

windws_snp_matrix_clean %>% 
  .[which(samps %in% samps_vec[grep("line", samps_vec,  invert = T)] ),
    which(snps_vec %in% guide_snp$SNP_id )
    ] ->
  in_matrix

D_cm <- dist(in_matrix)
tre_cm <- nj(D_cm)

tree_cm <- ggtree(tre_cm) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

is_tip <- tre_cm$edge[,2] <= length(tre_cm$tip.label)
ordered_tips <- tre_cm$edge[is_tip, 2]
tre_cm$tip.label[ordered_tips]  -> tree_order

joint_figure_polarized_missingdat_clean %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(samp_id %in%  samps_vec[grep("line", samps_vec,  invert = T)],
         SNP_id %in%  snps_vec ) %>% 
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

