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
### FST analysis

load("../2.Temporal_Spatial_structure/Year_to_year_object.Rdata")

Out_comp_vector_samepops %>%
  filter(bin_date == "2.Overwinter",
         month1 == month2,
         month1 %in% 7:11,
         pop1 == "Charlottesville") %>%
  ggplot(aes(
    x=as.factor(month1),
    y=FST
  )) + geom_boxplot(outlier.shape = NA) +
  geom_jitter() ->
  plot_box_fst_month

ggsave(plot_box_fst_month,
       file = "plot_box_fst_month.pdf")

Out_comp_vector_samepops %>%
  filter(bin_date == "2.Overwinter",
         month1 == month2,
         month1 %in% 7:11,
         pop1 == "Charlottesville") %>%
  group_by(month1) %>%
  summarise(FST_sd = sd(FST)) ->
  variance_of_FST
  
names(variance_of_FST)[1] = "month_col"

variance_of_FST %>%
  ggplot(aes(
    x=as.factor(month_col),
    y=FST_sd
  )) + geom_bar(stat = "identity")  ->
  plot_box_fst_month_var

ggsave(plot_box_fst_month_var,
       file = "plot_box_fst_month_var.pdf")


##########
### plot month
#####

## merge fst variance with pca table



######
######
pca_table_df %>%
  filter(year >= 2016) %>%
  group_by(month_col, chr) %>%
  summarise(sum_dim1 = sd(Dim.1),
            sum_dim2 = sd(Dim.2),
            sum_dim3 = sd(Dim.3)) %>% 
  left_join(variance_of_FST) ->
  var_per_month

var_per_month %>%
  filter(month_col %in% 7:11) %>% 
  melt(id = c("month_col", "chr", "FST_sd")) %>% 
  ggplot(aes(
    x=FST_sd,
    y=value,
    shape = chr,
    color = month_col
  )) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", 
              color = "black",
              se = F) +
  xlab(expression(Var(italic(F)[ST]))) +
  ylab("Var(PC projection)") +
  facet_grid(chr~variable) ->
  var_month

ggsave(var_month,
       file = "var_month.pdf")

var_per_month %>%
  filter(year > 2016,
         month_col %in% 7:11) ->
  data_in

lm(sum_dim2 ~ month_col + I(month_col^2), data = data_in[which(data_in$chr == "3R"),] ) ->
  lm_sq_month
summary(lm_sq_month)




