#### Enrrichemtn analysis

args = commandArgs(trailingOnly=TRUE)
start=args[1]
end=args[2]

setwd("/scratch/yey2sn/Overwintering_ms/4.GML_plots")

################
## load packages
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

## import gds
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")
## import metdata from glm + ld analyses
#glm_ld_outs <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers.txt")
glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"

load(glm.file)

glm.out %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) -> 
  glm.outliers

glm.outliers %>%
  .$SNP_id %>%
  unique() %>%
  data.frame(SNP_id=.) %>% 
  separate(SNP_id, remove = F, into = c("chr", "pos"), sep = "_") ->
  glm.outliers.unique
  

## create snpo guide file
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"))

snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>% 
  filter(SNP_id %in% glm.outliers$SNP_id ) ->
  snps.dt.out

## import the glm outlier + ld object
seqSetFilter(genofile, variant.id=snps.dt.out$variant.id)
subset_geno = genofile

####
snps.dt.out %>%
  .[start:end,] ->
  snps.dt.out.sb


#####
annotation_set <- foreach(i=1:dim(snps.dt.out.sb)[1],
                          .combine = "rbind")%do%{ ## open DO loop ## MAYBE CHANGE TO DOPAR?
                            
                            print(i)
                            subset_geno = genofile
                            
                            seqSetFilter(subset_geno, variant.id=snps.dt.out.sb$variant.id[i])
                            
                            tmp_annot <-  data.frame(annot = seqGetData(subset_geno, "annotation/info/ANN")[[2]])
                            
                            #extract annotation and parition into categories
                            tmp_annot %>%
                              separate(annot,
                                       into = c(
                                         "Allele",
                                         "Annotation",
                                         "Putative_impact",
                                         "Gene_Name",
                                         "Gene_ID",
                                         "Feature_type",
                                         "Feature_ID",
                                         "Transcript_biotype",
                                         "Rank",
                                         "HGVS.c",
                                         "HGVS.p",
                                         "cDNA_position",
                                         "CDS_position",
                                         "Protein_position",
                                         "Distance_to_feature",
                                         "Notes"),
                                       sep = "\\|",
                                       convert = T
                              ) -> tmp_annot_sep
                            
                            ## include annotaions as well as position and chromosome
                            obj_to_rturn <- data.frame(chr = snps.dt.out.sb$chr[i],
                                                       pos = snps.dt.out.sb$pos[i],
                                                       tmp_annot_sep
                            )
                            
                            return(obj_to_rturn) ##<---- return this object to the output object
                            
                          } ## close do loop


#####
#####

### Process annotations using a list decomposition format
### I will use a foreach loop to extract the annotation of each outlier SNP
annotation_set <- foreach(i=1:dim(snps.dt.out.sb)[1],
                          .combine = "rbind")%do%{ ## open DO loop ## MAYBE CHANGE TO DOPAR?
                            
                            print(i)
                            subset_geno = genofile
                            
                            seqSetFilter(subset_geno, variant.id=snps.dt.out.sb$variant.id[i])
                            
                            tmp_annot <-  data.frame(annot = seqGetData(subset_geno, "annotation/info/ANN")[[2]])
                            
                            #extract annotation and parition into categories
                            tmp_annot %>%
                              separate(annot,
                                       into = c(
                                         "Allele",
                                         "Annotation",
                                         "Putative_impact",
                                         "Gene_Name",
                                         "Gene_ID",
                                         "Feature_type",
                                         "Feature_ID",
                                         "Transcript_biotype",
                                         "Rank",
                                         "HGVS.c",
                                         "HGVS.p",
                                         "cDNA_position",
                                         "CDS_position",
                                         "Protein_position",
                                         "Distance_to_feature",
                                         "Notes"),
                                       sep = "\\|",
                                       convert = T
                              ) -> tmp_annot_sep
                            
                            ## include annotaions as well as position and chromosome
                            obj_to_rturn <- data.frame(chr = snps.dt.out.sb$chr[i],
                                                       pos = snps.dt.out.sb$pos[i],
                                                       tmp_annot_sep
                            )
                            
                            return(obj_to_rturn) ##<---- return this object to the output object
                            
                          } ## close do loop
#fix a boolean misinterpretation
annotation_set$Allele[which(annotation_set$Allele == "TRUE")] = "T"

### add priority scoring
### 
annotation_set %<>%
  mutate(priority_scoring = case_when(Annotation == "3_prime_UTR_variant" ~ 40,
                                      Annotation == "3_prime_UTR_premature_start_codon_gain_variant" ~ 45,
                                      Annotation == "5_prime_UTR_variant" ~ 40,
                                      Annotation == "downstream_gene_variant" ~ 0,
                                      Annotation == "intergenic_region" ~ 0,
                                      Annotation == "intron_variant" ~ 20,
                                      Annotation == "missense_variant" ~ 80,
                                      Annotation == "missense_variant&splice_region_variant" ~ 100,
                                      Annotation == "non_coding_transcript_exon_variant" ~ 15,
                                      Annotation == "splice_region_variant" ~ 70,
                                      Annotation == "splice_region_variant&intron_variant" ~ 70,
                                      Annotation == "synonymous_variant" ~ 50,
                                      Annotation == "upstream_gene_variant" ~ 0
  ))

annotation_set$priority_scoring[is.na(annotation_set$priority_scoring)] = 0

annotation_set %>% 
  group_by(chr, pos) %>%
  slice_max(priority_scoring, with_ties = F) ->
  annotation_set_prioritized


save(annotation_set_prioritized,
     file = paste("./enrrich_out/enrrichment",start,end,"set.Rdata", sep = ".")
     )
