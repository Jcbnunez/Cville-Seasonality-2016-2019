### Run Analysis of PCA
### 
rm(list = ls())

#load data
# This R object was premade in script 1 ==> 1.Import_GDStoR.r
data_in="/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA.object.all.Rdata"
inmeta="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata"

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
library(data.table)

### load data
load(data_in)
load(inmeta)

samps = samps_EFFCOV

samps$city = gsub("Charlotttesville","Charlottesville", samps$city)
samps$city = gsub("Odesa","Odessa", samps$city )
samps$city[grep("Yesiloz", samps$city )] = "Yesiloz"
samps$city[grep("Chornobyl", samps$city )] = "Chernobyl"
samps$city[grep("Kyiv", samps$city )] = "Kyiv"

#ectracts variace
PCA_object$eig[1,2] -> D1VE
PCA_object$eig[2,2] -> D2VE
PCA_object$eig[3,2] -> D3VE
PCA_object$eig[4,2] -> D4VE

#Plot PCA
## Plot object
PCA_object$ind$coord[,c(1:20)] %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>% 
  left_join(., samps ) -> PCA_table

PCA_table$collectionDate = as.Date(PCA_table$collectionDate, format = "%m/%d/%Y")
####
#### Make the centroid of each population
label.df_PCA <- PCA_table %>% 
  group_by(city) %>% 
  summarize(Dim.1 = mean(Dim.1), 
            Dim.2 = mean(Dim.2),
            Dim.3 = mean(Dim.3)) 

label.df_PCA.time <- PCA_table %>% 
  group_by(city, year) %>% 
  summarize(Dim.1 = mean(Dim.1), 
            Dim.2 = mean(Dim.2),
            Dim.3 = mean(Dim.3),
            Dim.4 = mean(Dim.4)) 


##run correlations wit time
PCA_table$city %>% unique -> cities_to_select

corr_list = list()
for(i in 1:length(cities_to_select)){
  print(i)
  
  pc_proj <- PCA_table$Dim.1[which(PCA_table$city == cities_to_select[i])]
  time <- PCA_table$year[which(PCA_table$city == cities_to_select[i])]
  
  tmp = cor.test(pc_proj, time)
  
  tmp_df = data.frame(p_val = tmp$p.value, 
                      cor = tmp$estimate,
                      city = cities_to_select[i],
                      sig_sym = ifelse(tmp$p.value < 0.05, "S", "NS") )
  
  corr_list[[i]] = tmp_df
  
}

corr_list_df = do.call(rbind, corr_list)

corr_list_df %<>%
  mutate(R2 = cor^2)

#colorset = c(brewer.pal(n = 12, name = "Set3"),
#             "darkgoldenrod1",
#             "brown2",
#             "hotpink4",
#             "olivedrab4"
#)

### time projection plot
label.df_PCA.time %>%
  group_by(city) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(mean_y = mean(year)) %>%
  mutate(endpt = ifelse(year > mean_y, "max", "min")) %>%
  melt(id=c("city","year","mean_y", "endpt")) %>% 
  dcast(city ~ variable+endpt ) %>%
  left_join(corr_list_df) ->
  label.df_PCA.time.min.max

  ggplot() + 
    geom_point(data  = PCA_table,
               shape = 21,
               aes(x=Dim.1, y=Dim.2, 
                   fill = city,
                   size = 3,
                   alpha = 0.8))  + 
    geom_text_repel(data = label.df_PCA, size = 3 ,  
                    box.padding = 0.5, max.overlaps = Inf,
                    alpha = 1, 
                    aes(x=Dim.1, 
                        y=Dim.2, 
                        label = city))  + 
    geom_segment(data = label.df_PCA.time.min.max,
                 size = 0.9,
                 aes(x = Dim.1_min, 
                     y = Dim.2_min, 
                     xend = Dim.1_max, 
                     yend = Dim.2_max
                     ),
                 arrow = arrow(length = unit(0.25, "cm"))) +
    geom_text_repel(data = label.df_PCA.time.min.max, size = 2 ,  
                    box.padding = 0.5, max.overlaps = Inf,
                    alpha = 1, 
                    aes(x=Dim.1_min, 
                        y=Dim.2_min, 
                        label = round(p_val, 3)))  + 
    theme_classic() +
    xlab(paste("PC 1 (", round(D1VE, 2), "% VE)" , sep = "")) +
    ylab(paste("PC 2 (", round(D2VE, 2), "% VE)", sep = "")) +
    scale_color_manual(values = c("red","blue")) +
    theme(legend.position = "none")->
  time_projection

ggsave(time_projection, 
       file = "time_projection.pdf",
       width = 4,
       height = 4)

### plot PCA_12 with time
##### RUN LM MODELS

  lm_inner_list = list()
  for(i in 1:20){
    
    lm(PCA_table[,i] ~ lat + long  + MeanEC,
       data = PCA_table) %>% 
      summary() -> tmp_1
    
    tmp_1$coefficients %>% 
      as.data.frame() %>%
      .[-1,] %>%
      mutate(PC = i,
             variable = rownames(.),
             Samples = "all_pops",
             r_sq = tmp_1$r.squared) -> tmp
    
    lm_inner_list[[i]] = tmp
    
  } #i

lm_dat_df = do.call(rbind, lm_inner_list)

lm_dat_df %<>% 
  mutate(P_bon = p.adjust(`Pr(>|t|)`, method = "bonferroni"))

lm_dat_df %>%
  filter(P_bon < 0.01) %>%
  select(PC, P_bon, variable)

