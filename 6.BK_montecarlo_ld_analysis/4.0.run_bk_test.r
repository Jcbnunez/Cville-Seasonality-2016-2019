## Do BK simulations --- R portion

args = commandArgs(trailingOnly=TRUE)
run_model=args[1]
global_iterator=args[2]


#load bk metadata
load("./bk_ld_frequencies.Rdata")
#bk_ld_frequencies

#load temp metadata
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"

# load DEST metadata
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata")
samps <- samps_EFFCOV
samps$collectionDate = as.Date(samps$collectionDate, format = "%m/%d/%Y")

#load libraries
library(tidyverse)
library(magrittr)
library(reshape2)
library(lme4)

# parition data
bk_ld_frequencies$snp_1 %>%
  unique() %>%
  .[complete.cases(.)] ->
  guide_focal_snps

### load in allele frequency data
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata")

o %>%
  t() %>%
  as.data.frame() %>%
  mutate(SNP_id = paste(rownames(.), "SNP", sep = "_")) ->
  o_t

outer_list = list()
for(i in 1:length(guide_focal_snps)){ #### open i
  
  bk_ld_frequencies %>%
    filter(snp_1 == guide_focal_snps[i]) ->
    master_obj
  
  inner_list = list()
  for(j in 1:dim(master_obj)[1]){ #### open j for master object

    o_t %>%
      filter(SNP_id == master_obj$snp_2[j]) %>% 
      select(!SNP_id) %>%
      t %>% 
      as.data.frame() %>%
      mutate(sampleId = rownames(.),
             snp_1 = guide_focal_snps[i],
             snp_2 = paste(colnames(.), "SNP", sep = "_") )  %>% 
      left_join(master_obj) %>%
      left_join(samps) %>%
      left_join(weather.ave) -> obj_for_glm
    
      names(obj_for_glm)[1] = "AF"
      obj_for_glm$AF = as.numeric(obj_for_glm$AF)
     
       obj_for_glm %<>%
        mutate(expected_p = (AF*present)+((1-AF)*absent)) 
      
      obj_for_glm$binom_samp = NA
      for(k in 1:dim(obj_for_glm)[1]){ #open sampling loop for binom
        
        obj_for_glm$binom_samp[k] = rbinom(1, 
                                    round(obj_for_glm$MeanEC[k],0), 
                                    obj_for_glm$expected_p[k] )
        } # close binom sampling loop
      
      obj_for_glm %<>%
      mutate(montecarlo_p = binom_samp/MeanEC) %>%
      mutate(af_montecarlo_nEff=round(montecarlo_p*MeanEC)/MeanEC)
    
      obj_for_glm %<>% 
        mutate(af_real_nEff=round(AF*MeanEC)/MeanEC)

      if(run_model == "real"){ #### Open real test
      model_temp <-  glm(af_real_nEff~I(aveTemp/10),
                         data=obj_for_glm,
                         family=binomial(),
                         weights=MeanEC)
      
      #estimate Nagelkerke-pseudo-R2
      real_output = data.frame(focal_snp = unique(obj_for_glm$snp_1),
                               affected_snp = unique(obj_for_glm$snp_2),
                               test_type = unique(obj_for_glm$test_type),
                               pseudo_r2 = (1 - (model_temp$deviance/model_temp$null.deviance)),
                               type = "real")
      
      inner_list[[j]] = real_output
      }
      
      if(run_model == "simulation"){ #### Open real test
      model_montecarlo <-  glm(af_montecarlo_nEff~I(aveTemp/10),
                             data=obj_for_glm,
                             family=binomial(),
                             weights=MeanEC)
    
      montecarlo_out <- data.frame(focal_snp = unique(obj_for_glm$snp_1),
                              affected_snp = unique(obj_for_glm$snp_2),
                              test_type = unique(obj_for_glm$test_type),
                              pseudo_r2 = (1 - (model_montecarlo$deviance/model_montecarlo$null.deviance)),
                              type = "montecarlo")
      
      inner_list[[j]] = montecarlo_out
  } # close simulation
  
    outer_list[[i]]  = do.call(rbind, inner_list )
  }# close j for master object
  
} #### close i

results_df = do.call(rbind, outer_list )

system("mkdir ./bk_out_sims")

save(results_df,
     file = paste("./bk_out_sims/",
                  "iteration.",
                  global_iterator,
                  ".output.Rdata", 
                  sep = ""))
