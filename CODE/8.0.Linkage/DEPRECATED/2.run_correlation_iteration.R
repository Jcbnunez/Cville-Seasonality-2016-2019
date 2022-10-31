##### ---> Correlated allele frequencies, pt 2, prepping for the matrix
##### 

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)

##
##

args = commandArgs(trailingOnly=TRUE)
iterator = as.numeric(args[1])

#load the data
load("./data_for_covariance_test.Rdata")

lenght_of_guide = dim(unique_combs_sep)[1]

guides = seq(from = 1, to = lenght_of_guide, by = 700)

length(guides) ## array = 1 - 8713

unique_combs_w = unique_combs_sep[guides[iterator]:guides[iterator+1],]  

###
out_list = list()
for(i in 1: 700){
  
    unique_combs_w$V1[i] -> focal_snp
    unique_combs_w$V2[i] -> test_snp
  
  dat_glm_p1 %>%
    .[,which(colnames(.) == focal_snp)] %>%
    melt(value.name = "AF_focal") %>% 
    mutate(sampleId = rownames(.))  -> focal_tmp
  
  dat_glm_p1 %>%
    .[,which(colnames(.) == test_snp)] %>%
    melt(value.name = "AF_test") %>% 
    mutate(sampleId = rownames(.))  -> test_tmp
  
  left_join(focal_tmp, test_tmp) -> joint_tmp
  
  cor.test(joint_tmp$AF_focal, joint_tmp$AF_test) -> tmp_out
  
  tmp_out <- data.frame(
    snp1 = focal_snp,
    snp2 = test_snp,
    chr1 = unique_combs_w$chr1[i],
    chr2 = unique_combs_w$chr2[i],
    pos1 = unique_combs_w$pos1[i],
    pos2 = unique_combs_w$pos2[i],
    bp_dist = unique_combs_w$bp_dist[i],
    correlation = tmp_out$estimate,
    p.value = tmp_out$p.value)
  
  out_list[[i]] =   tmp_out
}

out_df = do.call(rbind, out_list)

save(out_df,
     file = paste("corrOut/correlation_analysis",iterator, "Rdata", sep = "." ) )
