### make PCA from Cville chromosomes 
### 

rm(list = ls())
# Load packages

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
#install_github('tavareshugo/windowscanr')
library(windowscanr)

######
###### Import R objects

objects <- c(
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2R.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3L.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3R.ECfiltered.Rdata"
)

pca_table_list = list()
for(i in 1:length(objects)){
  
load(objects[i])
o %>% colnames() %>% .[1] -> lead_snp

lead_snp %>% data.frame(headsnp = .) %>% 
  separate(headsnp , into = c("CHR","POS")) -> lead_snp_guide

print(lead_snp_guide$CHR)

filtered_samps_for_analysis %>%
  filter(city == "Charlottesville",
         MeanEC > 30) %>% 
  .$sampleId -> select_samples

o %>%
  as.data.frame() %>% 
  filter(rownames(.) %in%  select_samples) %>% 
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  PCA_object

PCA_object$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.),
         chr = lead_snp_guide$CHR) %>%  
  left_join(., filtered_samps_for_analysis ) -> PCA_table

pca_table_list[[i]] = PCA_table
}

pca_table_df = do.call(rbind, pca_table_list)

### add some extra time metadata

pca_table_df %<>%
mutate(month_col = month(as.Date(collectionDate, 
                                 format = c("%m/%d/%Y"))))
  

#### general plot
#### Make the arrows

  ggplot() + 
  geom_point(shape =21, 
             size = 3,
             alpha = 0.9,
             data = pca_table_df,
             aes(
               x=Dim.1,
               y=Dim.2,
               fill = year)
             ) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 2016) +
  theme_bw() +
    ylab("PC2") +
    xlab("PC1") +
  facet_wrap(~chr) ->
  cville_chr_pca

ggsave(cville_chr_pca,
       file = "cville_chr_pca.pdf",
       width = 5,
       height = 4)


#### Part 2 run regressions against various predictors

pca_table_df$chr %>% unique -> chrs_vec

regressions_list1 = list()
regressions_list2 = list()
regressions_list3 = list()

for(i in 1:length(chrs_vec)){
  
  print(chrs_vec[i])
  
lm(Dim.1 ~ year + 
           MeanEC + 
           Month + 
           `In(2L)t`  +  
           `In(2R)Ns`  + 
           `In(3L)P` +
           `In(3R)C` +   
           `In(3R)Mo`+ 
           `In(3R)Payne`, 
   data = pca_table_df[which(pca_table_df$chr == chrs_vec[i]),] ) -> 
  lm_out

  summary(lm_out) -> lm_summ
  
  lm_summ$coefficients %>%
    as.data.frame() %>%
    mutate(variable = rownames(.),
           chr = chrs_vec[i],
           PC = 1,
           signif = ifelse(.$`Pr(>|t|)` < 0.05, "S", "NS")) ->
    lm_sum_fin

  regressions_list1[[i]] = lm_sum_fin
  
  #######

  lm(Dim.2 ~ year + 
       MeanEC + 
       Month + 
       `In(2L)t`  +  
       `In(2R)Ns`  + 
       `In(3L)P` +
       `In(3R)C` +   
       `In(3R)Mo`+ 
       `In(3R)Payne`, 
     data = pca_table_df[which(pca_table_df$chr == chrs_vec[i]),] ) -> 
    lm_out
  
  summary(lm_out) -> lm_summ
  
  lm_summ$coefficients %>%
    as.data.frame() %>%
    mutate(variable = rownames(.),
           chr = chrs_vec[i],
           PC = 2,
           signif = ifelse(.$`Pr(>|t|)` < 0.05, "S", "NS")) ->
    lm_sum_fin
  
  regressions_list2[[i]] = lm_sum_fin
  
  #####
  
  lm(Dim.3 ~ year + 
       MeanEC + 
       Month + 
       `In(2L)t`  +  
       `In(2R)Ns`  + 
       `In(3L)P` +
       `In(3R)C` +   
       `In(3R)Mo`+ 
       `In(3R)Payne`, 
     data = pca_table_df[which(pca_table_df$chr == chrs_vec[i]),] ) -> 
    lm_out
  
  summary(lm_out) -> lm_summ
  
  lm_summ$coefficients %>%
    as.data.frame() %>%
    mutate(variable = rownames(.),
           chr = chrs_vec[i],
           PC = 3,
           signif = ifelse(.$`Pr(>|t|)` < 0.05, "S", "NS")) ->
    lm_sum_fin
  
  regressions_list3[[i]] = lm_sum_fin
  
}

regressions_df <-
rbind(do.call(rbind, regressions_list1),
      do.call(rbind, regressions_list2),
      do.call(rbind, regressions_list3))

write.table(regressions_df, file = "./regression_cville_pcs.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


regressions_df %>% 
  filter(signif == "S",
         variable != "(Intercept)") %>% 
  select(chr, PC, variable, `Pr(>|t|)`) %>%
  dcast(PC+chr~variable)


##### Month Analysis

##########
### plot specific panels
#####

pca_table_df %>%
  filter(chr == "2L") %>%
  ggplot(aes(
    x=`In(2L)t`,
    y=Dim.1
  )) + 
  geom_point(size = 3) +
  geom_smooth(method = "lm") ->
  inv_vs_pc1

ggsave(inv_vs_pc1,
       file = "inv_vs_pc1.pdf",
       width = 4,
       height = 4)

## Plot variance pannel

pca_table_df %>%
  filter(chr != "2L") %>%
  ggplot(aes(
    x=year,
    y=Dim.2
  )) + 
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  facet_wrap(~chr, ncol =2) ->
  year_vs_pc2

ggsave(year_vs_pc2,
       file = "year_vs_pc2.pdf",
       width = 4,
       height = 4)

### projection panel

pca_table_df %>%
  ggplot(aes(
    x=as.factor(month_col),
    y=abs(Dim.1)
  )) +
  geom_boxplot() +
  geom_point() ->
  month_proejctions

ggsave(month_proejctions,
       file = "month_proejctions.pdf")

pca_table_df %>%
  ggplot(aes(
    x=as.factor(month_col),
    y=abs(Dim.2)
  )) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~chr)->
  month_proejction2

ggsave(month_proejction2,
       file = "month_proejction2.pdf")

