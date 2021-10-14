# Extract temperature loci from Temperature gml
# 

#load libraries
library(tidyverse)
library(magrittr)
library(data.table)

# add link to data
obj <- "/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata"
#load the object
load(obj)

#load annotations
load("/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata")

## extarct outlier SNPs
glm.out %>%
  filter(mod=="aveTemp",
         rnp.clean<0.01,
         chr=="2L") ->
  glm.out_p1

## add annotations
glm.out_p1 %>%
  left_join(annotation_dmel_sp, by = c("chr","pos")) ->
  glm.out_p1_annot

## add SNP annotation
glm.out_p1_annot %<>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))

# add ranking
glm.out_p1_annot %<>% 
  .[order(.$rnp.clean, decreasing = FALSE),] %>%
  mutate(final_glm_rank = 1:dim(.)[1])

## save_final_table
write.table(glm.out_p1_annot,
            file = "./Cville_glm_p01_hits_annot.txt",
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)

## save_temporal_table
write.table(glm.out_p1_annot$SNP_id[which(glm.out_p1_annot$final_glm_rank %in% 1:100)],
            file = "./temperature_snps_ids_top100.txt",
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)

# Load linkage SNPs for guide
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")

write.table(final_in2Lt_markers,
            file = "./final_in2Lt_markers.txt",
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)




