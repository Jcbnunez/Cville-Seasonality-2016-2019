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


#setwd("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/")

### load haplotype info

#load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata")
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
load("/scratch/yey2sn/Overwintering_ms/7.LD/merged.ld.Rdata")
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

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_9.6" ),
             mid = c(3.1, 4.7, 5.1, 6.1, 9.6)
             ) %>%
  mutate(start = mid-0.2 ,
          end  = mid+0.2  )

final.windows.pos

             #start = c(3000026, 4650065, 5050026, 6100321, 9500026),
             #end = c(3150026, 4799922,  5250026, 6224905, 9650026)) 


file.mod <- "/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"
glm.out <- get(load(file.mod))

glm.out %<>%
  mutate(fdr.score = p.adjust(p_lrt, method =  "fdr"))

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
head(final_in2Lt_markers)


glm.out %>%
  filter(perm == 0) %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(chr == "2L",
         fdr.score < 0.1) %>% 
  mutate(win = case_when(
    pos/1e6 > 3.1-0.2 & pos/1e6  < 3.1+0.2 ~ "win_3.1",
    pos/1e6  > 4.7-0.2 & pos/1e6  < 4.7+0.2 ~ "win_4.7",
    pos/1e6  > 5.1-0.2 & pos/1e6  < 5.1+0.2 ~ "win_5.1",
    pos/1e6  > 6.1-0.2 & pos/1e6  < 6.1+0.2 ~ "win_6.1",
    pos/1e6  > 9.6-0.2 & pos/1e6  < 9.6+0.2 ~ "win_9.6",
    pos < 3e6 & SNP_id %in% final_in2Lt_markers ~ "left",
    pos > 11e6 & SNP_id %in% final_in2Lt_markers ~ "right")) -> 
  glm.outliers.2L.wins

glm.outliers.2L.wins %>%
  filter(win %in% c("win_3.1","win_4.7","win_5.1","win_6.1","win_9.6"))  -> 
  glm.outliers.2L.wins.flt

rejected_pos = c("2L_3035162_SNP", "2L_3079654_SNP")
#good ones ->"2L_3056539_SNP", "2L_3056539_SNP"
  
glm.outliers.2L.wins.flt %>%
  filter(!SNP_id %in% rejected_pos  ) %>%
  group_by(win) %>%
  slice_min(fdr.score) %>% as.data.frame() %>% .$SNP_id ->
  top_snps

ld_df %>%
  filter(BP_A != BP_B) %>%
  #filter(SNP_A %in% top_snps | SNP_B %in% top_snps ) %>%
  filter(R2 > 0.60) %>% 
  #filter(BP_diff < 4e5) %>% 
  #group_by(pair_id) %>%
  #mutate(anchor = min(BP_A, BP_B)) %>% 
  #mutate(pair = case_when(anchor == BP_A ~ BP_B,
  #                       anchor != BP_A ~ BP_A)) %>% 
  mutate(win_B = 
  case_when(
    BP_B/1e6 > 3.1-0.2 & BP_B/1e6 < 3.1+0.2 ~ "win_3.1",
    BP_B/1e6 > 4.7-0.2 & BP_B/1e6 < 4.7+0.2 ~ "win_4.7",
    BP_B/1e6 > 5.1-0.2 & BP_B/1e6 < 5.1+0.2 ~ "win_5.1",
    BP_B/1e6 > 6.1-0.2 & BP_B/1e6 < 6.1+0.2 ~ "win_6.1",
    BP_B/1e6 > 9.6-0.2 & BP_B/1e6 < 9.6+0.2 ~ "win_9.6",
  )) %>%
  mutate(win_A = 
             case_when(
               BP_A/1e6 > 3.1-0.2 & BP_A/1e6 < 3.1+0.2 ~ "win_3.1",
               BP_A/1e6 > 4.7-0.2 & BP_A/1e6 < 4.7+0.2 ~ "win_4.7",
               BP_A/1e6 > 5.1-0.2 & BP_A/1e6 < 5.1+0.2 ~ "win_5.1",
               BP_A/1e6 > 6.1-0.2 & BP_A/1e6 < 6.1+0.2 ~ "win_6.1",
               BP_A/1e6 > 9.6-0.2 & BP_A/1e6 < 9.6+0.2 ~ "win_9.6",
             )) %>% 
  filter(win_A == win_B) ->
  #mutate(pos = as.numeric(BP_B) ) %>%
  #left_join(dplyr::select(sim_polarity, pos, tidy_annot)) ->
  haplo_tag_snps_r2
  
haplo_tag_snps_r2 %>% 
  group_by(win_B, SNP_B) %>% 
  summarize(N= n()) %>%
  group_by(win_B) %>%
  #slice_max(N, n = 10) %>%
  mutate(chr = "2L", SNP_id = SNP_B)-> haplotag_candidates

left_join(haplotag_candidates, glm.outliers.2L.wins.flt, by = c("chr", "SNP_id")) %>%
  .[complete.cases(.$pos),] %>% 
  group_by(win_B) %>%
  filter(fdr.score < 0.1) %>%
  slice_min(fdr.score, n = 1) %>% 
  as.data.frame() ->
  Anchor_SNPs

Anchor_SNPs

ld_df %>%
  filter(BP_A != BP_B) %>%
  filter(SNP_A %in% Anchor_SNPs$SNP_id | SNP_B %in% Anchor_SNPs$SNP_id ) %>%
  filter(R2 > 0.55) %>% 
  #filter(BP_diff < 4e5) %>% 
  #group_by(pair_id) %>%
  #mutate(anchor = min(BP_A, BP_B)) %>% 
  #mutate(pair = case_when(anchor == BP_A ~ BP_B,
  #                       anchor != BP_A ~ BP_A)) %>% 
  mutate(win_B = 
           case_when(
             BP_B/1e6 > 3.1-0.2 & BP_B/1e6 < 3.1+0.2 ~ "win_3.1",
             BP_B/1e6 > 4.7-0.2 & BP_B/1e6 < 4.7+0.2 ~ "win_4.7",
             BP_B/1e6 > 5.1-0.2 & BP_B/1e6 < 5.1+0.2 ~ "win_5.1",
             BP_B/1e6 > 6.1-0.2 & BP_B/1e6 < 6.1+0.2 ~ "win_6.1",
             BP_B/1e6 > 9.6-0.2 & BP_B/1e6 < 9.6+0.2 ~ "win_9.6",
           )) %>%
  mutate(win_A = 
           case_when(
             BP_A/1e6 > 3.1-0.2 & BP_A/1e6 < 3.1+0.2 ~ "win_3.1",
             BP_A/1e6 > 4.7-0.2 & BP_A/1e6 < 4.7+0.2 ~ "win_4.7",
             BP_A/1e6 > 5.1-0.2 & BP_A/1e6 < 5.1+0.2 ~ "win_5.1",
             BP_A/1e6 > 6.1-0.2 & BP_A/1e6 < 6.1+0.2 ~ "win_6.1",
             BP_A/1e6 > 9.6-0.2 & BP_A/1e6 < 9.6+0.2 ~ "win_9.6",
           )) %>% 
  filter(win_A == win_B) ->
  #mutate(pos = as.numeric(BP_B) ) %>%
  #left_join(dplyr::select(sim_polarity, pos, tidy_annot)) ->
  haplo_tag_snps_r2.final

save(haplo_tag_snps_r2.final, file = "haplo_tag_snps_r2.final.Rdata")
load("./haplo_tag_snps_r2.final.Rdata")


haplo_tag_snps_r2 %>%
  group_by(SNP_A, win_B) %>%
  summarize(Mean_R2 = mean(R2),
            Npairs = n()) -> haplotag_summaries

haplotag_summaries %>%
  filter(SNP_A %in% Anchor_SNPs$SNP_B ) ->
  #group_by(win_B) %>% 
  #filter(Npairs > 10) %>%
  #slice_max(Npairs, with_ties =F ) -> 
  Anchor_snps_of_haplotags

Anchor_snps_of_haplotags


glm.outliers.2L.wins %>%
  filter(SNP_id %in% Anchor_snps_of_haplotags$SNP_A)


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

