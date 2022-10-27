## This code is used to delineate the 2Lt inversion using the DGRP
## 

#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(RVenn)

#Metadata
invst <- fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
names(invst)[1] = "DGRP_Line"
invst$line_name = gsub("DGRP" , "line" , invst$DGRP_Line)

# Load the DGRP vcf
DGRP <- read.vcfR(
  "/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz")
DGRPgl <- vcfR2genlight(DGRP) 

#transform to genlight object
#check individuals
#tab(DGRPgl) %>%
#  rownames()

#Extract data for imputation
tab(DGRPgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  all_dgrp_data


#Separate data based on inversion set
standard  <- invst$line_name[which(invst$In_2L_t == "ST")]
inversion <- invst$line_name[which(invst$In_2L_t == "INV")]
heterokaryotype <- invst$line_name[which(invst$In_2L_t == "INV/ST")]

#Make standard set
all_dgrp_data %>%
  .[,which(colnames(.) %in% standard)] ->
  standard_lines

all_dgrp_data %>%
  .[,which(colnames(.) %in% inversion)] ->
  inverted_lines

all_dgrp_data %>%
  .[,which(colnames(.) %in% heterokaryotype)] ->
  heterokaryotype_lines

######## Define functions to study sets

filter_drgp_sets = function(x, miss_thresh){
  
  count_NA = function(x) {
    x  %>%   
      .[is.na(.)] %>%
      length() ->
      numNA
        return(numNA/length(x))
  }
  
  apply(x,
        2, 
        FUN=count_NA
  ) -> x_NA_lines
  
  x_keep_samps <- names(x_NA_lines[x_NA_lines <= miss_thresh])
  
  x %>%
    .[,which(colnames(.) %in%  x_keep_samps)] -> 
    x_lines_samps_filt
  
  #count NA per SNPs
  apply(x_lines_samps_filt,
        1, 
        FUN=count_NA
  ) -> x_NA_SNPs
  
  x_keep_snps <- names(x_NA_SNPs[x_NA_SNPs <= miss_thresh])
  
  x_lines_samps_filt %>%
    .[which(rownames(.) %in%  x_keep_snps),] -> 
    x_lines_samps_SNPs_filt
  
  x_lines_samps_SNPs_filt_NAimp = na.aggregate(x_lines_samps_SNPs_filt,
                                               FUN= median)
  
  return(x_lines_samps_SNPs_filt_NAimp)
  
}

######## ######## ######## ######## ######## 
######## Evaluate the standard_lines
standard_lines_flt = filter_drgp_sets(standard_lines, 1)

inverted_lines_flt = filter_drgp_sets(inverted_lines, 1)

heterokaryotype_lines_flt = filter_drgp_sets(heterokaryotype_lines, 1)

standard_lines_flt %>% dim
inverted_lines_flt %>% dim
heterokaryotype_lines_flt %>% dim

### standarize the datasets
standard_lines_flt %<>% 
  as.data.frame() 
  names(standard_lines_flt) = gsub("^", "STD_", names(standard_lines_flt))
  standard_lines_flt %<>% mutate(SNP_id = rownames(.))
  
inverted_lines_flt %<>% 
  as.data.frame()
  names(inverted_lines_flt) = gsub("^", "INV_", names(inverted_lines_flt))
  inverted_lines_flt %<>% mutate(SNP_id = rownames(.))
  
heterokaryotype_lines_flt %<>% 
  as.data.frame() 
  names(heterokaryotype_lines_flt) = gsub("^", "HET_", names(heterokaryotype_lines_flt))
  heterokaryotype_lines_flt %<>% mutate(SNP_id = rownames(.))
  
  ### merge the datasets
  ### Make metadta
  obj = list(
    std = rownames(standard_lines_flt),
    inverted_lines_flt = rownames(inverted_lines_flt),
    het = rownames(heterokaryotype_lines_flt)
  )

  venn_snps = Venn(obj)  
  overlap(venn_snps) -> shared_snps
  
  ### Merge by shared snps 
  standard_lines_flt %>%
    .[which(rownames(.) %in%  shared_snps), -which(names(.) == "SNP_id") ]->
    standard_lines_flt_shared

  inverted_lines_flt %>%
    .[which(rownames(.) %in%  shared_snps), -which(names(.) == "SNP_id") ]->
    inverted_lines_flt_shared

  heterokaryotype_lines_flt %>%
    .[which(rownames(.) %in%  shared_snps), -which(names(.) == "SNP_id") ]->
    heterokaryotype_lines_flt_shared

  
  cbind(standard_lines_flt_shared,
        inverted_lines_flt_shared,
        heterokaryotype_lines_flt_shared
        ) ->
    imputated_data_DGRP
  
#Save object
save(imputated_data_DGRP,
     file = "./imputated_data_DGRP_dat.2L.Rdata")
