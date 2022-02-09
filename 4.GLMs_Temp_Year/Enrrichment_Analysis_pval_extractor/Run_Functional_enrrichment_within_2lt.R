library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(doMC)
library(lubridate)

registerDoMC(8)

### This code incorporates the new method for enrrichent estimation based on fitting FETs
### internally for each type of dataset. It uses a soft match-control by restricting the
### unverse "N" of observartions to inversions/chromosomes sets. As opposed to the whole
### genome. This is done to ameliorate the inherent ascertainment bias of "large denominators"
### vs "small numerators"

args = commandArgs(trailingOnly=TRUE)

### Import some user given parameters
pop = args[1]
print(paste("now loading", pop))
#pop <- "VA_ch"
p_tresh = args[2]
#p_tresh = 0.05
print(paste("now processing p <", p_tresh))

#One could, also turn the model into an iterator to compare the T and Y+T models
#for now just Y+T
model = "aveTemp+year_factor"



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

## Clear up the annotation categories
## I do this to simplify the analysis
print("now cleaning annotation categories")
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
## I know that iteration 99 got bugged and never ran. This is ok, I will not iterate over 99
iterations = c(0,1:98,100)

#output lists
##Local_enrrichment_out = list()

###########################
### RUN THE FOREACH LOOP  #
###########################

Local_enrrichment_out <- foreach(i=  1:length(iterations) )%dopar%{
#for(i in 1:length(iterations)){ ### open i
  
  ### Part 1 build a metadata object with the genomic annotations
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
    filter(chr == "2L") %>%
    #filter(invName == "2Lt") %>%
    filter(!is.na(rnp.clean)) %>%
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
  
  #### Add annotatiobs
    glm.out_f %>%
    left_join(annotation_dmel_sp) ->
    glm.out_f_annotations
  
  glm.out_f_annotations %<>%
    .[complete.cases(.$Consequence),] %>%
    mutate(Consequence = case_when(Consequence %in% 
                                     c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence,CRM_pos, sep = "|"),
                                   !Consequence %in% 
                                     c("intergenic_region","intron_variant") ~ 
                                     paste(Consequence))) 
  ##
  glm.out_f_annotations %<>%
    filter(!Consequence %in% c("start","stop","NC_exon")) %>%
    mutate(Tidy_Consequence = case_when(
      Consequence == "3_prime" ~ "3p",
      Consequence == "5_prime" ~ "5p",
      Consequence == "intergenic_region|CRM" ~ "Inter+CRM",
      Consequence == "intergenic_region|non.CRM" ~ "Inter",
      Consequence == "intron_variant|CRM" ~ "Intron+CRM",
      Consequence == "intron_variant|non.CRM" ~ "Intron",
      Consequence == "missense_variant" ~ "Nonsyn",
      Consequence == "synonymous_variant" ~ "Syn",
      Consequence == "Splice_var" ~ "Splice",
      Consequence == "start" ~ "Start",
      Consequence == "stop" ~ "Stop",
      Consequence == "NC_exon" ~ "NC_ex"
    ))
  
  ### Make table of GLM oiutliers
  glm.out_f_annotations %>%
    filter(gene_type == "protein_coding") %>%
    group_by(Tidy_Consequence, Symbol) %>%
    summarize(All_SNPs = n()) -> 
    all_snps
  
  ### Make table of GLM oiutliers
  glm.out_f_annotations %>%
    filter(gene_type == "protein_coding") %>%
    filter(rnp.clean < p_tresh) %>% 
    group_by(Tidy_Consequence, Symbol) %>%
    summarize(Outliers = n()) -> 
    outlier_snps
  
  ## Merge tables
  full_join(all_snps, outlier_snps) -> 
    snp_table

  snp_table$Outliers[is.na(snp_table$Outliers)] = 0
  
  snp_table %>%
    mutate(proportion = Outliers/All_SNPs) %>% 
    mutate(chr = "2L",
           p_tresh = p_tresh,
           pop = pop,
           Analysis = "All",
           category = category,
           ith = iterations[i],
           ) ->
    all_snps
  
  
  #####
  ### Keep only functional loci
  
  analyses_types = list(
    NonSyn = c("Nonsyn"),
    Syn = c("Syn"),
    Coding = c("Syn", "Nonsyn"),
    Exon_Intron = c("Syn", "Nonsyn", "Intron", "Intron+CRM", "Splice", "3p", "5p"),
    Trascribed = c("Syn", "Nonsyn", "3p", "5p")
  )
  
  functional_out = list()
  
  for(k in  1:length(analyses_types) ){
    
    glm.out_f_annotations %>%
      filter(Tidy_Consequence %in% analyses_types[[k]] ) %>%
      group_by(Tidy_Consequence, Symbol) %>%
      summarize(All_SNPs = n()) -> 
      all_snps_tmp
    
    ### Make table of GLM oiutliers
    glm.out_f_annotations %>%
      filter(Tidy_Consequence %in% analyses_types[[k]] ) %>%
      filter(rnp.clean < p_tresh) %>% 
      group_by(Tidy_Consequence, Symbol) %>%
      summarize(Outliers = n()) -> 
      outlier_all_snps
    
    ## Merge tables
    full_join(all_snps_tmp, outlier_all_snps) -> 
      snp_table_tmp
    
    snp_table_tmp$Outliers[is.na(snp_table_tmp$Outliers)] = 0
    
    snp_table_tmp %>%
      mutate(proportion = Outliers/All_SNPs) %>% 
      mutate(chr = "2L",
             p_tresh = p_tresh,
             pop = pop,
             Analysis = names(analyses_types)[k],
             category = category,
             ith = iterations[i],
      ) ->
      tmp_out
   
    functional_out[[k]] = tmp_out
     
  } # close k
  
  
  functional_out_df = do.call(rbind, functional_out)
  
  return(rbind(all_snps, functional_out_df))
  
  #######
  
} ### Close i!!!

rbindlist(Local_enrrichment_out) -> Local_enrrichment_df

save(Local_enrrichment_df,
     file = paste(pop,p_tresh, "2L.gene.enrrich_count.Rdata", sep ="."))

