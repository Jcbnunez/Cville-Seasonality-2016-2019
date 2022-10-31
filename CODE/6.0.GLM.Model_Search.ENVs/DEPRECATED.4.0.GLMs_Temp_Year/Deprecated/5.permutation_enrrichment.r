#### 5. enrrichment of observed and permiutation data
####

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

###

#load annonations
load("/project/berglandlab/thermal_glm_dest_ANNOTATIONs/GLM_loci_annonations.Rdata")
#---> this loads annots_id

p_tresh=0.05

###
### Generate list of GLMs
glm_va_files <-
system("ls /project/berglandlab/thermal_glm_dest/processedGLM/ | grep 'VA' ", intern = T)

###
glms_out_merged = list()

#glm.out.VA_ch_42.Rdata
for(i in 1:length(glm_va_files)){
file <-
paste("/project/berglandlab/thermal_glm_dest/processedGLM/", 
      glm_va_files[i], 
      sep = "" )
# load file
message(paste(glm_va_files[i],"; i=" ,i))
load(file)
#---> glm.out
glm.out %>%
  filter(p.lrt < p_tresh) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) -> 
  glm.ids
### merge
left_join(glm.ids, annots_id) ->
  glm.ids_annot
glms_out_merged[[i]] = glm.ids_annot
} ## 

glms_out_merged_df = do.call(rbind, glms_out_merged)

#save(glms_out_merged_df, file = "enrrichment_analysis_w_perm.Rdata" )
load("/project/berglandlab/thermal_glm_dest_ANNOTATIONs/GLM_loci_annonations.Rdata")
#---> this loads annots_id

load("./enrrichment_analysis_w_perm.Rdata")
#--> glms_out_merged_df

###
###
###

total_snps = dim(annots_id)[1]
  
glms_out_merged_df %>% 
  filter(mod  == "aveTemp+year_factor") %>%
  group_by(perm, chr, invName, Annotation) %>%
  summarize(N = n()) %>% 
  #mutate(Tidy_Consequence = case_when(
  #  Annotation == "3_prime_UTR_variant" ~ "3p",
  #  Annotation == "5_prime_UTR_variant" ~ "5p",
  #  Annotation == "downstream_gene_variant" ~ "Inter",
  #  Annotation == "intron_variant" ~ "Intron",
  #  Annotation == "missense_variant" ~ "NS",
  #  Annotation == "missense_variant&splice_region_variant" ~ "NS",
  #  Annotation == "splice_acceptor_variant&intron_variant" ~ "Splice",
  #  Annotation == "splice_donor_variant&intron_variant" ~ "Splice",
  #  Annotation == "splice_region_variant" ~ "Splice",
  #  Annotation == "splice_region_variant&intron_variant" ~ "Splice",
  #  Annotation == "splice_region_variant&stop_retained_variant" ~ "Splice",
  #  Annotation == "splice_region_variant&synonymous_variant" ~ "Splice",
  #  Annotation == "start_lost" ~ "change_start_stop",
  #  Annotation == "stop_gained" ~ "change_start_stop",
  #  Annotation == "stop_lost" ~ "change_start_stop",
  #  Annotation == "stop_lost&splice_region_variant " ~ "change_start_stop",
  #  Annotation == "stop_retained_variant" ~ "change_start_stop",
  #  Annotation == "synonymous_variant" ~ "Syn",
  #  Annotation == "upstream_gene_variant" ~ "Inter",
  #  Annotation == "intergenic_region" ~ "Inter"
  #)) %>%
  mutate(Freq = N/total_snps) ->
  glms_out_merged_df_annot_tidy

glms_out_merged_df_annot_tidy %>%
  filter(Annotation %in% 
           c("3_prime_UTR_variant",
             "5_prime_UTR_variant",
             "intron_variant",
             "missense_variant",
             "synonymous_variant"
             ),
         chr != "X"
           ) -> dat_in

dat_in$Annotation = gsub("_variant", "", dat_in$Annotation)
dat_in$Annotation = gsub("_", " ", dat_in$Annotation)
dat_in$Annotation = gsub("prime", "'", dat_in$Annotation)

###
###
###

ggplot() +
  geom_violin(data = filter(dat_in, perm != 0),
              aes(
                x=chr,
                y=Freq,
                fill=invName)) +
  geom_point(data = filter(dat_in, perm == 0),
              aes(
                x=chr,
                y=Freq,
                fill=invName),
             shape = 23,
             position= position_dodge(width=0.91),
             size = 2.3
             ) +
  facet_wrap(chr~Annotation, scales = "free") ->
  conseq_violin_point

ggsave(conseq_violin_point, w= 9, h = 7, file = "conseq_violin_point.pdf")


  
  
