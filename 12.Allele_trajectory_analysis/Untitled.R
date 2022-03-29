
library(rehh)
library(patchwork)
library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(data.table)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(car)
library(DescTools)
library(viridis)
registerDoMC(20)

pheno_out = readRDS("/project/berglandlab/Adam/gwas.grouped")
load("SNP_guide_haps_wins_metadata.Rdata")
#----> SNP_guide_metadata


##pos > 5155762 & pos < 5255762 ~ "win5",
##pos > 6255762 & pos < 6355762 ~ "win6",
##pos > 9505762 & pos < 9605762 ~ "win9",
##left_break = getData(chr="2L", start=2051609 , end=3096574)
##right_break = getData(chr="2L", start=11584064 , end=13204668)


snps_enrriched = foreach(i=1:dim(SNP_guide_metadata)[1], .combine = "rbind", .errorhandling = "remove")%do%{
  
  snp_anchor = SNP_guide_metadata$pos[i]
  message(snp_anchor)
  pheno_out %>% 
    filter(grepl(  snp_anchor  , snp.positions)) %>%
    mutate(snp_anchor = snp_anchor )
  
}

save(snps_enrriched, file = "snps_enrriched.fullDF.Rdata")

snps_enrriched %>% 
  arrange(snp_anchor) %>%
  group_by(snp_anchor, gwas.pheno, group) %>%
  slice_head() %>%
  group_by(pos=snp_anchor) %>%
  summarise(N =n()) -> sum_tab

left_join(sum_tab, SNP_guide_metadata ) %>% 
  ggplot(aes(
    x=as.factor(loci_id),
    y=N,
    #fill = tidy_annot
  )) + geom_bar(stat = "identity") +
  facet_grid(.~win, scales = "free", space = "free"
  ) +
  #ggtitle(analysis) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")  ->
  pheno_plot
ggsave(pheno_plot, file =  "pheno_plot.pdf",  h = 4, w = 10)

  


