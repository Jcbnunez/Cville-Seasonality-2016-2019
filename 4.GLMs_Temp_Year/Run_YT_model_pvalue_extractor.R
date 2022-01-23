library(tidyverse)
library(magrittr)
#library(ggbeeswarm)

args = commandArgs(trailingOnly=TRUE)

### Import some user given parameters
pop = args[1]
#pop <- "VA_ch"

### where files are in Rivanna
annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
root_path <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out."

####### -- > model Y+T
######## 
## set some generalities
## Set in the path for annonations
## This steps take a bit of time
load(annotation_file)

## First declare the iterator
p_tresh = 0.05
iterations = c(0,1:98,100)
model = "aveTemp+year_factor"

#output lists
count_list= list()
Errich_list= list()

for(i in 1:length(iterations)){ ### open i

  print(iterations[i])
  
  load(paste(root_path,
             pop,
             "_", 
             iterations[i], 
             ".Rdata", sep = ""))
  
  category = ifelse(iterations[i] == 0, "Obs", "Per")
  print(category)
  
  glm.out %>%
    filter(chr %in% c("2L","2R","3L", "3R")) ->
    glm.out_f
  
  glm.out_f %>%
    filter(mod == model) %>%
    summarise(N=n()) -> all_snps
  
  glm.out_f %>%
    filter(mod == model) %>%
    filter(p.lrt < p_tresh) ->
    filtered_glm.out_f
  
  #count outliers
  filtered_glm.out_f %>%
  summarise(Tresh=n()) -> tresh_snp
  
  ## this step creates a data frame with the global counts
  genome_counts <- data.frame(chr="genome", 
                   tresh_snp,
                   all_snps,
                   ith = iterations[i],
                   p_tresh=p_tresh,
                   category=category,
                   pop=pop,
                   analysis_type = "global")
  
  #Next I will geenrate the counts based on chrosome
    glm.out_f %>%
    filter(mod == model) %>%
    group_by(chr) %>%
    summarise(N=n()) -> all_snps_chr
  
    filtered_glm.out_f %>%
    group_by(chr) %>%
    summarise(Tresh=n()) -> tresh_snp_chr
  
  chromosome_counts <- data.frame(tresh_snp_chr,
             all_snps_chr[-1],
             ith = iterations[i],
             p_tresh=p_tresh,
             category=category,
             pop=pop,
             analysis_type = "global")
  ### count analysis
  count_list[[i]] = rbind(genome_counts, chromosome_counts)
  
 ### Enrrichment analysis
 ### Firt merge with annotation object
  filtered_glm.out_f %>%
    left_join(annotation_dmel_sp) ->
    filtered_glm.out_f_annot
  
  ## Now generate inverions list
  filtered_glm.out_f_annot %<>%
    mutate(inversion_pos = case_when( chr == "2L" & pos > 2225744 & pos < 13154180  ~ "inv",
                                      chr == "2L" & pos < 2225744 | pos > 13154180  ~ "non.inv",
                                      chr == "2R" & pos > 15391154 & pos < 20276334  ~ "inv",
                                      chr == "2R" & pos < 15391154 | pos > 20276334  ~ "non.inv",
                                      chr == "3R" & pos > 11750567 & pos < 26140370  ~ "inv",
                                      chr == "3R" & pos < 11750567 | pos > 26140370  ~ "non.inv",
                                      chr == "3R" & pos > 21406917 & pos < 29031297  ~ "inv",
                                      chr == "3R" & pos < 21406917 | pos > 29031297  ~ "non.inv",
                                      chr == "3R" & pos > 16432209 & pos < 24744010  ~ "inv",
                                      chr == "3R" & pos < 16432209 | pos > 24744010  ~ "non.inv",
                                      chr == "3L" & pos > 3173046 & pos < 16308841  ~ "inv",
                                      chr == "3L" & pos < 3173046 | pos > 16308841  ~ "non.inv",
                                      ))
  
  ## Now generate simplify annotations for enrichment analysis
  filtered_glm.out_f_annot$Consequence[grep("3_prime", filtered_glm.out_f_annot$Consequence)] = "3_prime"
  filtered_glm.out_f_annot$Consequence[grep("5_prime", filtered_glm.out_f_annot$Consequence)] = "5_prime"
  filtered_glm.out_f_annot$Consequence[grep("splice", filtered_glm.out_f_annot$Consequence)] = "Splice_var"
  filtered_glm.out_f_annot$Consequence[grep("stop", filtered_glm.out_f_annot$Consequence)] = "stop"
  filtered_glm.out_f_annot$Consequence[grep("downstream", filtered_glm.out_f_annot$Consequence)] = "intergenic_region"
  filtered_glm.out_f_annot$Consequence[grep("upstream", filtered_glm.out_f_annot$Consequence)] = "intergenic_region"
  filtered_glm.out_f_annot$Consequence[grep("initiator", filtered_glm.out_f_annot$Consequence)] = "start"
  filtered_glm.out_f_annot$Consequence[grep("start_lost", filtered_glm.out_f_annot$Consequence)] = "start"
  filtered_glm.out_f_annot$Consequence[grep("non_coding_transcript_exon_variant", filtered_glm.out_f_annot$Consequence)] = "NC_exon"
  
  filtered_glm.out_f_annot$Consequence %>% table
  
  ## Now generate count tables
  filtered_glm.out_f_annot %>%
    group_by(chr, inversion_pos) %>%
    summarize(N = n()) %>%
    mutate(p_tresh=p_tresh,
           category=category,
           pop=pop,
           analysis_type = "enrr_Inv") ->
    count_by_Inv
  
  ## Now generate count tables by conseq
  filtered_glm.out_f_annot %>%
    group_by(chr, Consequence, inversion_pos) %>%
    summarize(N = n()) %>%
    mutate(p_tresh=p_tresh,
           category=category,
           pop=pop,
           analysis_type = "enrr_consq") ->
    count_by_Consequence_Inv
  
  ## Now generate count tables by effect
  filtered_glm.out_f_annot %>%
    group_by(chr, Effect, inversion_pos) %>%
    summarize(N = n()) %>%
    mutate(p_tresh=p_tresh,
           category=category,
           pop=pop,
           analysis_type = "enrr_Eff") %>%
    .[complete.cases(.),] ->
    count_by_Consequence_Eff
  
  names(count_by_Consequence_Eff)[2] = "Consequence"
  
  ## Merge tables
  rbind(count_by_Inv, count_by_Consequence_Inv, count_by_Consequence_Eff) ->
    enrrichment_table
  
  Errich_list[[i]] = enrrichment_table
  
} ### close i

count_df = do.call(rbind, count_list)
enrrich_df = do.call(rbind, Errich_list)

save(count_df, enrrich_df,
     file = paste(pop, ".enrrich_count.Rdata", sep =""))

##out_df %>%
##  as.data.frame() %>%
##  ggplot(aes(
##    x=chr,
##    y=Tresh/N,
##    color=category,
##    size=category
##  )) +
##  geom_beeswarm() +
##  scale_size_manual(values = c(4,1.5))-> 
##  p_tresh_plot
##
##ggsave(p_tresh_plot, file = "p_tresh_plot.pdf")
##
