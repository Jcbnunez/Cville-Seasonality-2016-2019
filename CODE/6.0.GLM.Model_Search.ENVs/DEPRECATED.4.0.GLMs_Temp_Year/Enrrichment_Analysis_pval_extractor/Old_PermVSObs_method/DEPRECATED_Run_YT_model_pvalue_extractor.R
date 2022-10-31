library(tidyverse)
library(magrittr)
library(data.table)
#library(ggbeeswarm)

args = commandArgs(trailingOnly=TRUE)

### Import some user given parameters
pop = args[1]
print(paste("now loading", pop))
#pop <- "VA_ch"
p_tresh = args[2]
#p_tresh = 0.05
print(paste("now processing p <", p_tresh))

### where files are in Rivanna
annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
regulatory_map <- "/project/berglandlab/Dmel_genomic_resources/RegulatoryMap_database/remap2022_crm_macs2_dm6_v1_0.bed"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

root_path <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out."

####### -- > model Y+T
######## 
## set some generalities
## Set in the path for annonations
## This steps take a bit of time
print("now loading annotation file")
load(annotation_file)


annotation_dmel_sp$Consequence[grep("3_prime", annotation_dmel_sp$Consequence)] = "3_prime"
annotation_dmel_sp$Consequence[grep("5_prime", annotation_dmel_sp$Consequence)] = "5_prime"
annotation_dmel_sp$Consequence[grep("splice", annotation_dmel_sp$Consequence)] = "Splice_var"
annotation_dmel_sp$Consequence[grep("stop", annotation_dmel_sp$Consequence)] = "stop"
annotation_dmel_sp$Consequence[grep("downstream", annotation_dmel_sp$Consequence)] = "intergenic_region"
annotation_dmel_sp$Consequence[grep("upstream", annotation_dmel_sp$Consequence)] = "intergenic_region"
annotation_dmel_sp$Consequence[grep("initiator", annotation_dmel_sp$Consequence)] = "start"
annotation_dmel_sp$Consequence[grep("start_lost", annotation_dmel_sp$Consequence)] = "start"
annotation_dmel_sp$Consequence[grep("non_coding_transcript_exon_variant", annotation_dmel_sp$Consequence)] = "NC_exon"

print("now loading regulatory map file")
rmap<- fread(regulatory_map)

print("now loading inversion map file")
inv_map <- fread(inversion_map)
names(inv_map)[1:3] =  c("chr","start","end")
setkey(inv_map, chr, start, end)

######
rmap_simple = rmap[,1:3] 
names(rmap_simple) = c("chr","start","end")
rmap_simple %<>%
  mutate(regulatory = "CRM")

rmap_simple$chr = gsub("chr", "", rmap_simple$chr)
setkey(rmap_simple, chr, start, end)


## First declare the iterator
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
    filter(mod == model) %>%
    filter(chr %in% c("2L","2R","3L", "3R")) %>%
    mutate(start = pos, 
           end = pos)->
    glm.out_f
  
  ### Build overlaps
  foverlaps(glm.out_f[,c("chr","start","end")], inv_map, nomatch=NA) %>%
    mutate(inversion_pos = case_when( invName %in%  inv_map$invName ~ "inv",
                                      is.na(invName) ~ "non.inv")) %>%
    dplyr::select(chr, pos=i.start, inversion_pos) %>%
    mutate(SNP_id = paste(chr,pos, sep = "_")) %>% 
    distinct(SNP_id, .keep_all = TRUE) ->
    inversion_forverlaps
  
  foverlaps(glm.out_f[,c("chr","start","end")], rmap_simple, nomatch=NA) %>% 
    mutate(CRM_pos = case_when( regulatory == "CRM" ~ "CRM",
                                      is.na(regulatory) ~ "non.CRM")) %>%
    dplyr::select(chr, pos=i.start, CRM_pos) %>%
    mutate(SNP_id = paste(chr,pos, sep = "_")) %>% 
    distinct(SNP_id, .keep_all = TRUE) ->
    CRM_forverlaps
  
  ## sanity check
  print("overlap check -- must see two 'TRUE' statements")
  print(dim(glm.out_f)[1] == dim(CRM_forverlaps)[1])
  print(dim(glm.out_f)[1] == dim(inversion_forverlaps)[1])
  
  ####
  ####
  glm.out_f %<>%
  left_join(inversion_forverlaps[,c("chr","pos", "inversion_pos")] ) %>% 
  left_join(CRM_forverlaps[,c("chr","pos", "CRM_pos")] ) 
    
  ####
  ####
  glm.out_f %>%
    summarise(N=n()) -> all_snps
  
  glm.out_f %>%
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
  
  #Next I will generate the counts based on chromosome
    glm.out_f %>%
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
  
  ###Break by inversion status
  glm.out_f %>%
    group_by(chr, inversion_pos) %>%
    summarise(N=n()) -> all_snps_chr_inv
  
  filtered_glm.out_f %>%
    group_by(chr, inversion_pos) %>%
    summarise(Tresh=n()) -> tresh_snp_chr_inv
  
  chromosome_inv_counts <- data.frame(tresh_snp_chr_inv,
                                  all_snps_chr_inv[-c(1:2)],
                                  ith = iterations[i],
                                  p_tresh=p_tresh,
                                  category=category,
                                  pop=pop,
                                  analysis_type = "global_inversion")

  ###Break by CRM status
  glm.out_f %>%
    group_by(chr, CRM_pos) %>%
    summarise(N=n()) -> all_snps_chr_CRM
  
  filtered_glm.out_f %>%
    group_by(chr, CRM_pos) %>%
    summarise(Tresh=n()) -> tresh_snp_chr_CRM
  
  chromosome_CRM_counts <- data.frame(tresh_snp_chr_CRM,
                                      all_snps_chr_CRM[-c(1:2)],
                                      ith = iterations[i],
                                      p_tresh=p_tresh,
                                      category=category,
                                      pop=pop,
                                      analysis_type = "global_CRM")
  
  
  ### Joint INV/CRM
  glm.out_f %>%
    group_by(chr, CRM_pos, inversion_pos) %>%
    summarise(N=n()) -> all_snps_chr_CRM_INV
  
  filtered_glm.out_f %>%
    group_by(chr, CRM_pos, inversion_pos) %>%
    summarise(Tresh=n()) -> tresh_snp_chr_CRM_INV
  
  chromosome_CRM_INV_counts <- data.frame(tresh_snp_chr_CRM_INV,
                                      all_snps_chr_CRM_INV[-c(1:3)],
                                      ith = iterations[i],
                                      p_tresh=p_tresh,
                                      category=category,
                                      pop=pop,
                                      analysis_type = "global_CRM_INV")
  
  
  ### Joint INV/CRM/Functional Cats
  glm.out_f %>%
    left_join(annotation_dmel_sp) ->
    glm.out_f_annotations
  
  glm.out_f_annotations %<>%
    .[complete.cases(.$Consequence),] %>%
    mutate(Consequence = case_when(Consequence %in% c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence,CRM_pos, sep = "|"),
                                   !Consequence %in% c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence))) 
  
  #glm.out_f_annotations$Consequence %>% table
  
  glm.out_f_annotations %<>%
    filter(!Consequence %in% c("start","stop","NC_exon","Splice_var")) %>%
    mutate(Tidy_Consequence = case_when(
      Consequence == "3_prime" ~ "3p",
      Consequence == "5_prime" ~ "5p",
      Consequence == "intergenic_region|CRM" ~ "Inter+CRM",
      Consequence == "intergenic_region|non.CRM" ~ "Inter",
      Consequence == "intron_variant|CRM" ~ "Intron+CRM",
      Consequence == "intron_variant|non.CRM" ~ "Intron",
      Consequence == "missense_variant" ~ "Nonsyn",
      Consequence == "synonymous_variant" ~ "Syn"
    ))
    
  glm.out_f_annotations$Tidy_Name %>% table
  
  glm.out_f_annotations %>%
    group_by(chr, Tidy_Consequence, inversion_pos) %>%
    summarise(N=n()) -> all_snps_chr_CRM_INV_Annot
  
  glm.out_f_annotations %>%
    filter(p.lrt < p_tresh) %>% 
    group_by(chr, Tidy_Consequence, inversion_pos) %>%
    summarise(Tresh=n()) -> Tresh_chr_CRM_INV_Annot
  
  all_snps_chr_CRM_INV_Annot_counts <-  data.frame(Tresh_chr_CRM_INV_Annot,
                                                   all_snps_chr_CRM_INV_Annot[-c(1:3)],
                                          ith = iterations[i],
                                          p_tresh=p_tresh,
                                          category=category,
                                          pop=pop,
                                          analysis_type = "Annotation_enrichment")
  
  
  ### count analysis
  rbind_obj <- bind_rows(genome_counts, 
              chromosome_counts, 
              chromosome_inv_counts, 
              chromosome_CRM_counts, 
              chromosome_CRM_INV_counts,
              all_snps_chr_CRM_INV_Annot_counts) #%>%
    #mutate(analysis_subtype = paste(chr, inversion_pos, CRM_pos, sep = ""))
    #rbind_obj$analysis_subtype = gsub("NA","", rbind_obj$analysis_subtype)
  
    count_list[[i]] = rbind_obj
  
  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ### 
 ### Enrrichment analysis
 ### Firt merge with annotation object
  filtered_glm.out_f %>%
    left_join(annotation_dmel_sp) ->
    filtered_glm.out_f_annot
  
  ## Now generate simplify annotations for enrichment analysis
  
  filtered_glm.out_f_annot %<>%
    .[complete.cases(.$Consequence),] %>%
    mutate(Consequence = case_when(Consequence %in% c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence,CRM_pos, sep = "|"),
                                   !Consequence %in% c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence))) 
  
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
     file = paste(pop,p_tresh, "enrrich_count.Rdata", sep ="."))
