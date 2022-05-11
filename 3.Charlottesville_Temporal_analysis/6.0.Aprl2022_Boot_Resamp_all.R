
### Run PCA on DEST data
### 
rm(list = ls())

#load data
# This R object was premade in script 1 ==> 1.Import_GDStoR.r
data_in="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.mAF_Miss_Mean_Filt.ECfiltered.Rdata"
name="all"

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(data.table)
library(foreach)
library(scales)


### user params
args = commandArgs(trailingOnly=TRUE)

array_num = args[1]
sample_nsnps = args[2]
cov = args[3]

message(paste(array_num, sample_nsnps, cov, sep = " -- "))

# Examples
#> array_num=1
#> sample_nsnps=100

root="/scratch/yey2sn/Overwintering_ms/3.Cville_PCA/Aprl2022_boot_resamp"



#### TEST FUNCTION

test.corrs.pca = function(Pop,chr,type,jobid,sample_nsnps, cov, object){
  
  out.df = 
  data.frame(
    Pop=Pop,
    chr=chr,
    type=type,
    jobid=jobid,
    cov=cov,
    sample_nsnps = sample_nsnps,
    
    Year_1 = cor.test(~Dim.1 + year, object)$est,
    Year_2 = cor.test(~Dim.2 + year, object)$est,
    Year_3 = cor.test(~Dim.3 + year, object)$est,
    
    In2Lt_1 = cor.test(~Dim.1 + `In(2L)t`, object)$est,
    In2Lt_2 = cor.test(~Dim.2 + `In(2L)t`, object)$est,
    In2Lt_3 = cor.test(~Dim.3 + `In(2L)t`, object)$est,
    
    In2Rns_1 = cor.test(~Dim.1 + `In(2R)Ns`, object)$est,
    In2Rns_2 = cor.test(~Dim.2 + `In(2R)Ns`, object)$est,
    In2Rns_3 = cor.test(~Dim.3 + `In(2R)Ns`, object)$est,
    
    In3Lp_1 = cor.test(~Dim.1 + `In(3L)P`, object)$est,
    In3Lp_2 = cor.test(~Dim.2 + `In(3L)P`, object)$est,
    In3Lp_3 = cor.test(~Dim.3 + `In(3L)P`, object)$est,
    
    In3Rc_1 = cor.test(~Dim.1 + `In(3R)C`, object)$est,
    In3Rc_2 = cor.test(~Dim.2 + `In(3R)C`, object)$est,
    In3Rc_3 = cor.test(~Dim.3 + `In(3R)C`, object)$est,
    
    In3RMo_1 = cor.test(~Dim.1 + `In(3R)Mo`, object)$est,
    In3RMo_2 = cor.test(~Dim.2 + `In(3R)Mo`, object)$est,
    In3RMo_3 = cor.test(~Dim.3 + `In(3R)Mo`, object)$est,
    
    In3Rpayne_1 = cor.test(~Dim.1 + `In(3R)Payne`, object)$est,
    In3Rpayne_2 = cor.test(~Dim.2 + `In(3R)Payne`, object)$est,
    In3Rpayne_3 = cor.test(~Dim.3 + `In(3R)Payne`, object)$est,
    
    NC_1 = cor.test(~Dim.1 + MeanEC, object)$est,
    NC_2 = cor.test(~Dim.2 + MeanEC, object)$est,
    NC_3 = cor.test(~Dim.3 + MeanEC, object)$est)
  
  return(out.df)
}


### load data
load(data_in)
### Prepare chromsome samples
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata")
head(samps_EFFCOV)
### This object contains extra samps for analysis -- a more detailed metadata file
### bring number of inversions
#load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata")
#filtered_samps_for_analysis %>%
#  dplyr::select(sampleId,`In(2L)t`,   `In(2R)Ns`, `In(3L)P`,   `In(3R)C` , `In(3R)Mo`,  `In(3R)Payne`) ->
#  inversion_freqs
#save(inversion_freqs, file = "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/inversion_freqs.Rdata")
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/inversion_freqs.Rdata")


selected_pops.chr = list(
  Cha = "VA_ch",
  Bro = "DE_Bro",
  Mun = "DE_Mun",
  #Gim = "ES_Gim",
  Aka = "FI_Aka",
  li = "PA_li",
  Yes = "TR_Yes",
  #Kyi= "UA_Kyi",
  Ode = c("UA_Ode", "UA_od"),
  Cp = "WI_cp"
)

#SNP_sampler <- c(seq(from = 100, to = 1000, by = 100),
#                 seq(from = 1000, to = 20000, by = 1000))


chrs = c("2L","2R","3L", "3R")

### REAL CORRS
correlations.test = foreach(k=1:length(chrs)
                            , .combine = "rbind", .errorhandling = "remove")%:% ## open outer nested loop 
  foreach(i=1:length(selected_pops.chr), .combine = "rbind", .errorhandling = "remove")%do%{ ## open inner loop
    
    
    colnames(dat_for_Analysis) -> all_snps
    snps_chr = all_snps[grep(chrs[k], all_snps)]
    
    message(paste(names(selected_pops.chr)[i],chrs[k], sep = "-in-" ))
             
  if(selected_pops.chr[[i]] == "VA_ch"){
    samps_EFFCOV %>%
      filter(locality %in% selected_pops.chr[[i]],
             MeanEC > 30)%>% 
      .$sampleId -> selected_samps
    
    cov = 30
    
    message("VA filter applied")
  }
  
    if(selected_pops.chr[[i]] != "VA_ch"){
  samps_EFFCOV %>%
      filter(locality %in% selected_pops.chr[[i]],
             MeanEC > 28)%>% 
      .$sampleId -> selected_samps
    
      cov = 28
      
      message("DEST filter preserved")
    }
             message(selected_samps)
             
             ### each chr
             dat_for_Analysis %>%
               as.data.frame() %>% 
               filter(rownames(.) %in% selected_samps) %>% 
               .[,which(colnames(.) %in% snps_chr)] -> o
             
             o %>%
               .[,sample(dim(o)[2],sample_nsnps )] %>% 
               PCA(scale.unit = F, graph = F, ncp = 3) ->
               PCA_object
             
             PCA_object$ind$coord %>%
               as.data.frame() %>%
               mutate(sampleId = rownames(.),
                      chr=chrs[k],
                      analysis_set = names(selected_pops.chr)[i] ) %>%
               left_join(samps_EFFCOV) %>%
               left_join(inversion_freqs) -> tmp_obj
             
             ##tmp_obj %>%
             ##  ggplot(aes(
             ##    x=Dim.1,
             ##    y=Dim.2,
             ##    fill=year
             ##  ))  + 
             ##  scale_fill_gradientn(
             ##    colors=c("springgreen","cyan","blue","gold","red"),
             ##    values=rescale(c(2009,2010,2013,2016,2018))
             ##  ) +
             ##  geom_point(shape = 21) ->
             ##  pca_plot.tmp
             ##ggsave(pca_plot.tmp, file = "pca_plot.tmp.pdf")
  
             #Pop,chr,type,jobid,sample_nsnps, cov, object
             test.corrs.pca(Pop=names(selected_pops.chr)[i],
                            chr=chrs[k],
                            type="Real",
                            jobid=array_num,
                            sample_nsnps = sample_nsnps,
                            cov=cov,
                            tmp_obj) -> real_tmp
   
             ## Begin shuffle
             ## Begin shuffle
             ## Begin shuffle
             ## Begin shuffle
             ## Begin shuffle
             
             tmp_shuff = tmp_obj
             
             tmp_shuff$Dim.1 = tmp_shuff$Dim.1[shuffle(tmp_shuff$Dim.1)]
             tmp_shuff$Dim.2 = tmp_shuff$Dim.2[shuffle(tmp_shuff$Dim.2)]
             tmp_shuff$Dim.3 = tmp_shuff$Dim.3[shuffle(tmp_shuff$Dim.3)]
             
             test.corrs.pca(
                 Pop=names(selected_pops.chr)[i],
                 chr=chrs[k],
                 type="Permutation",
                 jobid=array_num,
                 sample_nsnps = sample_nsnps,
                 cov=cov,
                 tmp_shuff) -> perm_tmp
                 
             rbind(real_tmp, perm_tmp)
              
  } ## close inner loop

write.table(correlations.test, 
            file = paste(root,"/boot_resample_",array_num, "_" , sample_nsnps, "_", cov,".txt", sep =""),
            append = FALSE,
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = TRUE)
