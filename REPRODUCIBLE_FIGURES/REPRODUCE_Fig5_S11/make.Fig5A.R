### Reproduce Panel 5A
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
load("./haplowins_pt1.Rdata")
###
rownames(windws_snp_matrix_clean) = windws_snp_matrix_clean$samp_id

analyses_types = list(
  all=c("left",
        "win_3.1",
        "win_4.7",
        "win_5.1",
        "win_6.1",
        "win_6.6",
        "win_9.6",
        "right") #,
  #win3.1=c("win_3.1"),
  #win5.1=c("win_5.1"),
  #win9.6=c("win_9.6")
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
  
  tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_missingdat_clean + geom_tippoint(aes(color=pop))
  
  
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



