#### Concatenate 

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(foreach)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

setwd("/scratch/yey2sn/Overwintering_ms/4.GML_plots")

files <- system("ls ./enrrich_out", intern = T)

annots <- foreach(i = 1:length(files), .combine = "rbind")%do%{
  message(files[i])
  load(paste( "./enrrich_out" , files[i], sep = "/") )
  tmp <- annotation_set_prioritized
  return(tmp)
}

annots %>% 
  mutate(SNP_id = paste(chr, pos, sep = "_")) -> 
  annots_id

#save(annots_id, file = "/project/berglandlab/thermal_glm_dest_ANNOTATIONs/GLM_loci_annonations.Rdata")

annots_id %>% dim

### bing in the outliers
glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"

load(glm.file)

glm.out %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) -> 
  glm.ids

glm.ids %>%
  filter(mod == "aveTemp+year_factor") %>%
  dim

### merge
left_join(glm.ids, annots_id) ->
  glm.ids_annot

glm.ids_annot %>%
  filter(Transcript_biotype == "protein_coding") %>% 
  mutate(Tidy_Consequence = case_when(
    Annotation == "3_prime_UTR_variant" ~ "3p",
    Annotation == "5_prime_UTR_variant" ~ "5p",
    Annotation == "downstream_gene_variant" ~ "Inter",
    Annotation == "intron_variant" ~ "Intron",
    Annotation == "missense_variant" ~ "NS",
    Annotation == "missense_variant&splice_region_variant" ~ "NS",
    Annotation == "splice_acceptor_variant&intron_variant" ~ "Splice",
    Annotation == "splice_donor_variant&intron_variant" ~ "Splice",
    Annotation == "splice_region_variant" ~ "Splice",
    Annotation == "splice_region_variant&intron_variant" ~ "Splice",
    Annotation == "splice_region_variant&stop_retained_variant" ~ "Splice",
    Annotation == "splice_region_variant&synonymous_variant" ~ "Splice",
    Annotation == "start_lost" ~ "change_start_stop",
    Annotation == "stop_gained" ~ "change_start_stop",
    Annotation == "stop_lost" ~ "change_start_stop",
    Annotation == "stop_lost&splice_region_variant " ~ "change_start_stop",
    Annotation == "stop_retained_variant" ~ "change_start_stop",
    Annotation == "synonymous_variant" ~ "Syn",
    Annotation == "upstream_gene_variant" ~ "Inter"
  )) -> glm.ids_annot_tidy

glm.ids_annot_tidy %>%
  .[complete.cases(.$Tidy_Consequence),] %>%
  filter(mod == "year_factor",
         chr != "X") %>%
  group_by(Tidy_Consequence, chr, invName) %>%
  summarise(total_n_snps = n()) %>%
  as.data.frame() %>%
  dplyr::select(Tidy_Consequence, chr, invName, total_n_snps) ->
  total_snps

p_tresh=0.05

glm.ids_annot_tidy %>%
  .[complete.cases(.$Tidy_Consequence),] %>%
  filter(p.lrt < p_tresh,
         chr != "X") %>% 
  group_by(mod, Tidy_Consequence, chr, invName) %>%
  summarise(outliers = n()) %>% 
  as.data.frame() %>% 
  reshape2::dcast(Tidy_Consequence+chr+invName~mod) -> 
  outlier_snps

left_join(outlier_snps, total_snps) ->
  merged_table_counts

merged_table_counts %>%
  ggplot(aes(
    x=`aveTemp+year_factor`/total_n_snps,
    y=year_factor/total_n_snps,
    color =Tidy_Consequence,
    shape=chr
  )) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  #xlim(0,0.5) +
  #ylim(0,0.5) +
  facet_wrap(~invName) -> 
  props_counts

ggsave(props_counts, file = "props_counts.pdf")

#####
##ORs

OR_tab = foreach(i=1:dim(merged_table_counts)[1], .combine = "rbind")%do%{
  
  #A --> outlier feature
  #B --> feature not outlier
  #C --> outlier not feature
  #D --> not feature not outlier
  
  merged_table_counts[is.na(merged_table_counts)] = 0
  
  A=merged_table_counts$`aveTemp+year_factor`[i]
  B=merged_table_counts$total_n_snps[i]-merged_table_counts$`aveTemp+year_factor`[i]
  C=merged_table_counts %>% 
    .[-i,] %>%
    summarise(sum(`aveTemp+year_factor`)) %>%
    .[1,1]
  D=merged_table_counts %>% 
    .[-i,] %>%
    summarise(sum(total_n_snps)) %>%
    .[1,1]
  

  matrix <-
    matrix(c(A, B, C, D),
           nrow = 2)

  #matrix
  
  fisher.test(matrix) -> fet_tmp
  
  out <- data.table(
    A=A,
    B=B,
    C=C,
    D=D,
    feature=merged_table_counts$Tidy_Consequence[i],
    chr=merged_table_counts$chr[i],
    invName=merged_table_counts$invName[i],
    p.val=fet_tmp$p.value,
    OR=fet_tmp$estimate,
    OR_l=fet_tmp$conf.int[1],
    OR_h=fet_tmp$conf.int[2],
    model="temperature_OR"
  )
  
  return(out)
  
}


#OR_tab %>% 
#  as.data.frame() %>%
#  filter(feature != "change_start_stop",
#         feature != "Splice") %>%
#  ggplot(
#    aes(
#      x=chr,
#      y=log2(OR),
#      ymin=log2(OR_l),
#      ymax=log2(OR_h),
#      color = invName
#    )) +
#  geom_hline(yintercept = 0) +
#  geom_errorbar(width=0.1, position=position_dodge(width=0.5)) +
#  geom_point( position=position_dodge(width=0.5) ) +
#  ylim(-2,2) +
#  facet_wrap(~feature) ->
#  OR_plot
#
#ggsave(OR_plot, file = "OR_plot.pdf", w = 6, h = 3)

###


OR_time = foreach(i=1:dim(merged_table_counts)[1], .combine = "rbind")%do%{
  
  #A --> outlier feature
  #B --> feature not outlier
  #C --> outlier not feature
  #D --> not feature not outlier
  
  merged_table_counts[is.na(merged_table_counts)] = 0
  
  A=merged_table_counts$year_factor[i]
  B=merged_table_counts$total_n_snps[i]-merged_table_counts$year_factor[i]
  C=merged_table_counts %>% 
    .[-i,] %>%
    summarise(sum(year_factor)) %>%
    .[1,1]
  D=merged_table_counts %>% 
    .[-i,] %>%
    summarise(sum(total_n_snps)) %>%
    .[1,1]
  
  
  matrix <-
    matrix(c(A, B, C, D),
           nrow = 2)
  
  #matrix
  
  fisher.test(matrix) -> fet_tmp
  
  out <- data.table(
    A=A,
    B=B,
    C=C,
    D=D,
    feature=merged_table_counts$Tidy_Consequence[i],
    chr=merged_table_counts$chr[i],
    invName=merged_table_counts$invName[i],
    p.val=fet_tmp$p.value,
    OR=fet_tmp$estimate,
    OR_l=fet_tmp$conf.int[1],
    OR_h=fet_tmp$conf.int[2],
    model="time_OR"
  )
  
  return(out)
  
}

####

#OR_time %>% 
#  as.data.frame() %>%
#  filter(feature != "change_start_stop",
#         feature != "Splice") %>%
#  ggplot(
#    aes(
#      x=chr,
#      y=log2(OR),
#      ymin=log2(OR_l),
#      ymax=log2(OR_h),
#      color = invName
#    )) +
#  geom_hline(yintercept = 0) +
#  geom_errorbar(width=0.1, position=position_dodge(width=0.5)) +
#  ggtitle("Time Model") +
#  geom_point( position=position_dodge(width=0.5) ) +
#  ylim(-2,2) +
#  facet_wrap(~feature) ->
#  OR_plot_t
#
#ggsave(OR_plot_t, file = "OR_plot_time.pdf", w = 6, h = 3)

###

OR_tab_yt = foreach(i=1:dim(merged_table_counts)[1], .combine = "rbind")%do%{
  
  
  merged_table_counts[is.na(merged_table_counts)] = 0
  
  
  A=merged_table_counts$`aveTemp+year_factor`[i]
  B=merged_table_counts$total_n_snps[i]-merged_table_counts$`aveTemp+year_factor`[i]
  
  C=merged_table_counts$`year_factor`[i]
  D=merged_table_counts$total_n_snps[i]-merged_table_counts$`year_factor`[i]
  

  matrix <-
    matrix(c(A, B, C, D),
           nrow = 2)
  
  #matrix
  
  fisher.test(matrix) -> fet_tmp
  
  out <- data.table(
    A=A,
    B=B,
    C=C,
    D=D,
    feature=merged_table_counts$Tidy_Consequence[i],
    chr=merged_table_counts$chr[i],
    invName=merged_table_counts$invName[i],
    p.val=fet_tmp$p.value,
    OR=fet_tmp$estimate,
    OR_l=fet_tmp$conf.int[1],
    OR_h=fet_tmp$conf.int[2],
    model="Temp/Time_OR"
    
  )
  
  return(out)
  
}

###

rbind(
OR_time,
OR_tab,
OR_tab_yt) %>%
  as.data.frame() %>%
  filter(feature != "change_start_stop",
         feature != "Splice") %>%
  ggplot(
    aes(
      x=chr,
      y=log2(OR),
      ymin=log2(OR_l),
      ymax=log2(OR_h),
      color = invName
    )) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width=0.1, position=position_dodge(width=0.5)) +
  ggtitle("3 Models") +
  geom_point( position=position_dodge(width=0.5) ) +
  #ylim(-2.5,1.6) +
  facet_grid(model~feature) ->
  merged_ORs

ggsave(merged_ORs, file = "merged_ORs.pdf", w = 9, h = 4)

###

##
##