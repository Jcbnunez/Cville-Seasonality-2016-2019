## Do BK simulations --- R portion
## 

args = commandArgs(trailingOnly=TRUE)

iterator=as.numeric(args[1])

print(iterator)

#load bk metadata
load("./bk_ld_frequencies.Rdata")

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
bk_ld_frequencies$focal_snp %>%
  unique() %>%
  .[which(. != "")] ->
  guide_focal_snps


##### subset LD set to selected SNP

selected_focal_snp = guide_focal_snps[iterator]


### load in allele frequency data
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata")

o %>%
  t() %>%
  as.data.frame() %>%
  mutate(SNP_id = paste(rownames(.), "SNP", sep = "_")) ->
  o_t

## select SNP of interest
bk_ld_frequencies %>%
  filter(focal_snp == selected_focal_snp) ->
  master_obj

master_obj %<>%
  filter(affected_snp %in% o_t$SNP_id)

####### begin loop

real_r2 =list()
montecarlo_r2 =list()
for(i in 1:dim(master_obj)[1]){

  o_t %>%
    filter(SNP_id == master_obj$affected_snp[i]) %>% 
    select(!SNP_id) %>%
    t %>% 
    as.data.frame() %>%
    mutate(sampleId = rownames(.),
           focal_snp = selected_focal_snp,
           affected_snp = paste(colnames(.), "SNP", sep = "_") )  %>% 
    left_join(master_obj) %>%
    left_join(samps) %>%
    left_join(weather.ave) -> tmp1
  
  names(tmp1)[1] = "AF"
  tmp1$AF = as.numeric(tmp1$AF)
  
  tmp1 %<>%
    mutate(expected_p = (AF*present)+((1-AF)*absent)) 
  
  
  tmp1 %<>% 
    mutate(af_real_nEff=round(AF*MeanEC)/MeanEC)

  
  model_temp <-  glm(af_real_nEff~I(aveTemp/10),
                     data=tmp1,
                     family=binomial(),
                     weights=MeanEC)
  
  #estimate Nagelkerke-pseudo-R2
  tmp_2 = data.frame(focal_snp = unique(tmp1$focal_snp),
                     affected_snp = unique(tmp1$affected_snp),
                     pseudo_r2 = (1 - (model_temp$deviance/model_temp$null.deviance)),
                     type = "real"
                     )
  
  real_r2[[i]] = tmp_2
  
  #### begin montecarlo
  #montecarlo_r2 =list()
  
  loop_montecarlo_out = list()
  for(j in 1:100){
    
  tmp1$binom_samp = NA
  for(k in 1:dim(tmp1)[1]){
    
  tmp1$binom_samp[k] = rbinom(1, 
                                round(tmp1$MeanEC[k],0), 
                                tmp1$expected_p[k] )
  } # close binom sampling loop
  
  tmp1 %<>%
  mutate(montecarlo_p = binom_samp/MeanEC) %>%
  mutate(af_montecarlo_nEff=round(montecarlo_p*MeanEC)/MeanEC)
  
  model_montecarlo <-  glm(af_montecarlo_nEff~I(aveTemp/10),
                     data=tmp1,
                     family=binomial(),
                     weights=MeanEC)
  
  tmp_3 <- data.frame(focal_snp = unique(tmp1$focal_snp),
                     affected_snp = unique(tmp1$affected_snp),
                     pseudo_r2 = (1 - (model_montecarlo$deviance/model_montecarlo$null.deviance)),
                     type = "montecarlo")
  
  loop_montecarlo_out[[j]] = tmp_3
  
  } # close loop of bk sim
  
  montecarlo_r2[[i]] = do.call(rbind, loop_montecarlo_out)

  } # close larger loop
 

 
montecarlo_df = do.call(rbind, montecarlo_r2)
real_df = do.call(rbind, real_r2)

rbind(real_df, montecarlo_df) -> dat_in

system("mkdir ./montecarlo_out")

save(dat_in,
  file = paste("./montecarlo_out/",
           selected_focal_snp,
           ".Rdata", 
           sep = ""))

## prob of explain by hithchikers
dat_in$affected_snp %>% unique -> affected_snps
out_p_list = list()
for(z in 1:length(affected_snps)){
  
  tmp_p <-    
    dat_in %>%
    summarize(p_val = 
                mean(dat_in[which(dat_in$affected_snp == affected_snps[z] & dat_in$type == "montecarlo" ),]$pseudo_r2 >= dat_in[which(dat_in$affected_snp == affected_snps[z] & dat_in$type == "real" ),]$pseudo_r2 ))
  
  out_p_list[[z]] = data.frame(affected_snps=affected_snps[z], p_val= tmp_p)
}

out_p_df = do.call(rbind,out_p_list)

system("mkdir ./montecarlo_pvals_out")

save(out_p_df,
     file = paste("./montecarlo_pvals_out/",
                  selected_focal_snp,
                  ".pval.Rdata", 
                  sep = ""))


## plots
## 

##ggplot() +
##  geom_violin(
##    data = dat_in[which(dat_in$type == "montecarlo"),],
##    color = "steelblue",
##    aes(
##      x=as.factor(affected_snp),
##      y=pseudo_r2,
##      color = type)) +
##  geom_point(
##    data = dat_in[which(dat_in$type == "real"),],
##    fill = "red",
##    shape = 5,
##    aes(
##      x=as.factor(affected_snp),
##      y=pseudo_r2,
##      color = type)) +
##  coord_flip() +
##  ggtitle(paste("Focal to", unique(tmp1$focal_snp, sep = " ") ),
##          subtitle = paste(unique(master_obj$V1))) ->
##  test_f
##
##ggsave(test_f, file = "test_f.pdf")

#out_p_df %>%
#  ggplot(aes(
#    x=affected_snps.z.,
#    y=p_val
#  )) + 
#  geom_point() +
#  coord_flip() +
#  ggtitle("Prob. that Temp GLM is explained by LD to focal marker alone") +
#  ylim(0,1) ->
#  probz
#
#ggsave(probz, file = "probz.pdf")

