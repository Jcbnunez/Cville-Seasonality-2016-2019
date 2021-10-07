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

pca_table_df %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = year
  )) + 
  geom_point(shape =21, 
             size = 3,
             alpha = 0.9) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 2016) +
  theme_bw() +
  facet_wrap(~chr) ->
  cville_chr_pca

ggsave(cville_chr_pca,
       file = "cville_chr_pca.pdf",
       width = 5,
       height = 4)


### plot month

pca_table_df %>%
  filter(year >= 2016) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = month_col,
    shape = as.factor(year)
  )) + 
  geom_point() +
  scale_color_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 9) +
  facet_wrap(~chr) ->
  cville_chr_pca_m

ggsave(cville_chr_pca_m,
       file = "cville_chr_pca_m.pdf")
#####


pca_table_df %>%
  filter(year >= 2016) %>% 
  select( Dim.1, Dim.2, year, chr) %>%
  melt(id = c("chr","year")) %>% 
  ggplot(aes(
    x=year,
    y=value,
    color = variable,
  )) + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~chr) ->
  cville_chr_pca_melt

ggsave(cville_chr_pca_melt,
       file = "cville_chr_pca_melt.pdf")

#####
pca_table_df %>%
  filter(year >= 2016) %>%
  group_by(month_col, chr) %>%
  summarise(sum_dim1 = sd(Dim.1),
            sum_dim2 = sd(Dim.2)) ->
  var_per_month

var_per_month %>%
  melt(id = c("month_col", "chr")) %>%
  ggplot(aes(
    x=month_col,
    y=value,
    color = variable
  )) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~chr) ->
  var_month

ggsave(var_month,
       file = "var_month.pdf")

lm(sum_dim2 ~ month_col + I(month_col^2), data = var_per_month[which(var_per_month$chr == "3R"),] ) ->
  lm_sq_month
summary(lm_sq_month)







pca_table_df %<>%
  left_join(var_per_month)

pca_table_df %>%
  filter(year >= 2016) %>%
  ggplot(aes(
    y=abs(Dim.1),
    x=as.factor(month_col),
    col = chr,
  )) +
  #geom_point(size = 3) +
  geom_boxplot() +
  facet_wrap(~chr, scales = "free") ->
  var_dim_vs_month

ggsave(var_dim_vs_month,
       file = "var_dim_vs_month.pdf")



###

