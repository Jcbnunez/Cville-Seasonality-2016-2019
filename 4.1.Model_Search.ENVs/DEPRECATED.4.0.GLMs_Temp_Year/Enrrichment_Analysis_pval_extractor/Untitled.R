library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)

load("./VA_ch.0.01.2L.gene.enrrich_count.Rdata")

genes <- unique(Local_enrrichment_df$Symbol)

##dat_in <- Local_enrrichment_df 

annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
load(annotation_file)

### Make gene addresses
annotation_dmel_sp %>%
  group_by(chr, Symbol) %>%
  summarise(Median_pos = median(pos)) ->
  position_genes

### Generate Id list
prot_cod_ids <- unique(Local_enrrichment_df$Symbol)

Local_enrrichment_df %>%
  group_by(Symbol, category, ith, Analysis) %>%
  summarise(All_SNPs_tot = sum(All_SNPs),
            Outliers_tot = sum(Outliers)) %>% 
  mutate(proportion = Outliers_tot/All_SNPs_tot) ->
  Local_enrrichment_df_summarize

###
Local_enrrichment_df_summarize %>%   
  left_join(position_genes) -> 
  Local_enrrichment_df_summarize_pos

Local_enrrichment_df_summarize_pos %<>% 
  group_by(category, Symbol, Median_pos, Analysis, All_SNPs_tot) %>%
  summarise(Median = median(proportion),
            prop_l = quantile(proportion, 0.025),
            prop_h = quantile(proportion, 0.975),
            prop_sd = sd(proportion)
            ) 

Local_enrrichment_df_summarize_pos %>%
  filter(category == "Per") %>% 
  dplyr::select(Median_pos, category, Analysis, Symbol, top_9X=prop_h) ->
  top_9X

left_join(Local_enrrichment_df_summarize_pos, top_9X[,-2]) -> 
  Local_enrrichment_df_summarize_pos_top

#save(Local_enrrichment_df_summarize_pos_top, file = "/project/berglandlab/jcbnunez/Shared_w_Alan/Local_enrrichment_df_summarize_pos_top.Rdata")

############## Bring in Windows!!
output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")
ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)
#### load windows
#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
load(output_results_window)
###
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0,
                       list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
                       list(mod, chr.x=chr.x, locality, win.i, start, end)] %>%
  filter(locality == "VA_ch",
         mod=="aveTemp+year_factor")
win.out %<>% filter(locality == "VA_ch", mod=="aveTemp+year_factor")
############## Bring in Windows!! -0- close



dat_in = Local_enrrichment_df_summarize_pos_top %>%
  mutate(Cat_pval = 1-pbinom(Median*All_SNPs_tot, All_SNPs_tot, .01)) %>%
  filter(All_SNPs_tot >= 5)
  #.[complete.cases(.$Cat_pval),]




  ggplot() +
  #geom_point(
  #  data=dat_in[which(dat_in$category == "Per"),],
  #  aes(
  #  x=Median_pos,
  #  y=-1*log10(Cat_pval),
  #  #ymin = prop_l,
  #  #ymax = prop_h,
  #  fill = category,
  #  #color = category
  #), alpha = 0.3) +
    geom_point(
      data=dat_in[which(dat_in$category == "Obs"),],
      aes(
        x=Median_pos,
        y=-1*log10(Cat_pval),
        fill = category), 
      size = 1,
      shape = 23
  )  +
    geom_point(
      data=dat_in[which(dat_in$category == "Obs" &
                        dat_in$Median > dat_in$top_9X),],
      aes(
        x=Median_pos,
        y=-1*log10(Cat_pval)), 
      fill = "firebrick",
      size = 1.5,
      shape = 23
    )  +
    geom_point(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
               aes(x=(start/2 +end/2) , y=1*log10(rbinom.p)/10,
                   color=(rnp.pr)), size=.95) +
    scale_color_viridis(option="G", direction = -1) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~Analysis) +
    ylab(expression( paste(italic(P[binomial])) ) ) +
    xlab("bp") +
    ggtitle("Enrrichment of functional loci") +
    geom_vline(xintercept = 2225744) + 
    geom_vline(xintercept = 13154180) -> 
  p_tresh_plot_genes

ggsave(p_tresh_plot_genes, 
       file = "p_tresh_plot_genes.pdf",
       width = 9,
       height = 4)
 

####
#### NOTES FOR ALAN
## 1. Load object
load("/project/berglandlab/jcbnunez/Shared_w_Alan/Local_enrrichment_df_summarize_pos_top.Rdata")

#what is in this object?
names(Local_enrrichment_df_summarize_pos_top)

#category: Obs or Per, observed or permuted
#Symbol: Gene name
#Median_pos: Median position of the Gene in the GFF
#Analysis: All SNPs, all tranlated snos, Coding SNPs, Syn, NonSyn
#All_SNPs_tot: Number of total SNPs that fall in this category
#Median: This is the proportion of SNPs < P 5% / All SNPs. For the "Per" category this is the median over 100 perms. For "Obs" is the real value
#prop_l: output of prop_l = quantile(proportion, 0.01)
#prop_h: output of prop_l = quantile(proportion, 0.95)
#prop_sd: sd of permutations
#top_9X: Value of the top 95% of permutations this conects Obs and Per