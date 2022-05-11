### Run Analysis of PCA
### 
rm(list = ls())

#load data
data_in = "/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA_obj_ind_analysis.allnalyses.Rdata"

#data_in="/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA.object.all.Rdata"
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
library(scales)

### load data
load(data_in)
load(inmeta)

#samps = samps_EFFCOV

PCA_obj_ind_analysis$city = gsub("Charlotttesville","Charlottesville", PCA_obj_ind_analysis$city)
PCA_obj_ind_analysis$city = gsub("Odesa","Odessa", PCA_obj_ind_analysis$city )
PCA_obj_ind_analysis$city[grep("Yesiloz", PCA_obj_ind_analysis$city )] = "Yesiloz"
#PCA_obj_ind_analysis$city[grep("Chornobyl", PCA_obj_ind_analysis$city )] = "Chernobyl"
#PCA_obj_ind_analysis$city[grep("Kyiv", PCA_obj_ind_analysis$city )] = "Kyiv"

#ectracts variace
##PCA_object$eig[1,2] -> D1VE
##PCA_object$eig[2,2] -> D2VE
##PCA_object$eig[3,2] -> D3VE
##PCA_object$eig[4,2] -> D4VE

#Plot PCA
## Plot object
PCA_obj_ind_analysis -> PCA_table
#  filter(analysis_set == "all") 


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
  
  time <- PCA_table$year[which(PCA_table$city == cities_to_select[i])]
  
  
  pc_proj1 <- PCA_table$Dim.1[which(PCA_table$city == cities_to_select[i])]
  tmp1 = cor.test(pc_proj1, time)
  tmp_df1 = data.frame(p_val = tmp1$p.value, 
                      cor = tmp1$estimate,
                      PC = 1,
                      city = cities_to_select[i],
                      sig_sym = ifelse(tmp1$p.value < 0.05, "S", "NS") )
  
  pc_proj2 <- PCA_table$Dim.2[which(PCA_table$city == cities_to_select[i])]
  tmp2 = cor.test(pc_proj2, time)
  tmp_df2 = data.frame(p_val = tmp2$p.value, 
                       cor = tmp2$estimate,
                       PC = 2,
                       city = cities_to_select[i],
                       sig_sym = ifelse(tmp2$p.value < 0.05, "S", "NS") )
  
  pc_proj3 <- PCA_table$Dim.3[which(PCA_table$city == cities_to_select[i])]
  tmp3 = cor.test(pc_proj3, time)
  tmp_df3 = data.frame(p_val = tmp3$p.value, 
                       cor = tmp3$estimate,
                       PC = 3,
                       city = cities_to_select[i],
                       sig_sym = ifelse(tmp3$p.value < 0.05, "S", "NS") )
  
  
  corr_list[[i]] = rbind(tmp_df1, tmp_df2, tmp_df3)
  
}
corr_list_df = do.call(rbind, corr_list)


corr_list_df %<>%
  mutate(R2 = cor^2)

corr_list_df %>%
  group_by(PC, sig_sym) %>%
  summarise(N = n())
  

corr_list_df

write.table(corr_list_df, file = "corr_list_df.DEST.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



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

####
ggplot() + 
  geom_point(data  = PCA_table,
             shape = 21,
             aes(x=Dim.1, y=Dim.2, 
                 fill = year,
                 size = 2.5,
                 alpha = 0.8))  + 
  scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
  geom_text(data = label.df_PCA, size = 3 ,  
                  box.padding = 0.5, #max.overlaps = Inf,
                  alpha = 1, 
                  aes(x=Dim.1, y=Dim.2, 
                      label = city))  + 
  geom_path(data = label.df_PCA.time, 
            aes(x=Dim.1, 
                y=Dim.2,
                group = city), arrow = arrow(length=unit(0.10,"cm"),  type = "closed"),size = 0.7) +
  #geom_text_repel(data = label.df_PCA.time, size = 1.2 ,  
  #                box.padding = 0.5, max.overlaps = Inf,
  #                alpha = 1, 
  #                aes(x=Dim.1, y=Dim.2, 
  #                    label = year))  + 
  #geom_segment(data = label.df_PCA.time.min.max,
  #             size = 0.9,
  #             aes(x = Dim.1_min, 
  #                 y = Dim.3_min, 
  #                 xend = Dim.1_max, 
  #                 yend = Dim.3_max
  #             ),
  #             arrow = arrow(length = unit(0.25, "cm"))) +
  #geom_text_repel(data = label.df_PCA.time.min.max, size = 2 ,  
  #                box.padding = 0.5, max.overlaps = Inf,
  #                alpha = 1, 
#                aes(x=Dim.1_min, 
#                    y=Dim.2_min, 
#                    label = round(p_val, 3)))  + 
theme_classic() +
  #xlab(paste("PC 1 (", round(D1VE, 2), "% VE)" , sep = "")) +
  #ylab(paste("PC 2 (", round(D2VE, 2), "% VE)", sep = "")) +
  scale_color_manual(values = c("red","blue")) +
  theme(legend.position = "left")->
  time_projection12

ggsave(time_projection12, 
       file = "time_projection12.pdf",
       width = 5,
       height = 4)


#### PC 2 and 3

ggplot() + 
  geom_point(data  = PCA_table,
             shape = 21,
             aes(x=Dim.2, y=Dim.3, 
                 fill = year,
                 size = 2.5,
                 alpha = 0.8))  + 
  scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
    geom_text(data = label.df_PCA, size = 3 ,  
                  box.padding = 0.5, #, max.overlaps = Inf,
                  alpha = 1, 
                  aes(x=Dim.2, y=Dim.3, 
                      label = city))  + 
  geom_path(data = label.df_PCA.time, 
            aes(x=Dim.2, 
                y=Dim.3,
                group = city), arrow = arrow(length=unit(0.10,"cm"),  type = "closed"),size = 0.7) +
  #geom_text_repel(data = label.df_PCA.time, size = 1.2 ,  
  #                box.padding = 0.5, max.overlaps = Inf,
  #                alpha = 1, 
  #                aes(x=Dim.2, y=Dim.3, 
  #                    label = year))  + 
  #geom_segment(data = label.df_PCA.time.min.max,
  #             size = 0.9,
  #             aes(x = Dim.1_min, 
  #                 y = Dim.3_min, 
  #                 xend = Dim.1_max, 
  #                 yend = Dim.3_max
  #             ),
  #             arrow = arrow(length = unit(0.25, "cm"))) +
  #geom_text_repel(data = label.df_PCA.time.min.max, size = 2 ,  
  #                box.padding = 0.5, max.overlaps = Inf,
  #                alpha = 1, 
  #                aes(x=Dim.1_min, 
  #                    y=Dim.2_min, 
  #                    label = round(p_val, 3)))  + 
  theme_classic() +
  #xlab(paste("PC 1 (", round(D1VE, 2), "% VE)" , sep = "")) +
  #ylab(paste("PC 2 (", round(D2VE, 2), "% VE)", sep = "")) +
  scale_color_manual(values = c("red","blue")) +
  theme(legend.position = "left")->
  time_projection23

ggsave(time_projection23, 
       file = "time_projection23.pdf",
       width = 5,
       height = 4)


ggsave(time_projection12+time_projection23, 
       file = "time_projection12_23.pdf",
       width = 10,
       height = 4)



### Multiyear plot
PCA_obj_ind_analysis %>%
  filter(analysis_set != "all",
         analysis_set != "ch",) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=year,
    label=year
  )) +
  geom_text_repel(size = 3.5)+
  geom_point(size = 2, shape = 21) +
  scale_fill_gradient2(midpoint = 2013, low = "blue", high = "red", 
                       mid = "gold") +
  facet_wrap(~city) ->
  multipop_plot
ggsave(multipop_plot, file = "multipop_plot.pdf")


### plot PCA_12 with time
##### RUN LM MODELS

  lm_inner_list = list()
  for(i in 1:20){
    
    lm(PCA_table[,i] ~ year + lat + long  + MeanEC,
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
    
    # Eu
    lm(PCA_table[which(PCA_table$continent == "Europe"),i] ~ 
         year + lat + long  + MeanEC,
       data = PCA_table[which(PCA_table$continent == "Europe"),]) %>% 
      summary() -> tmp_2
    
    tmp_2$coefficients %>% 
      as.data.frame() %>%
      .[-1,] %>%
      mutate(PC = i,
             variable = rownames(.),
             Samples = "Europe",
             r_sq = tmp_2$r.squared) -> tmp_eu
    
    # N. Ame
    lm(PCA_table[which(PCA_table$continent == "NorthAmerica"),i] ~ 
         year + lat + long  + MeanEC,
       data = PCA_table[which(PCA_table$continent == "NorthAmerica"),]) %>% 
      summary() -> tmp_3
    
    tmp_3$coefficients %>% 
      as.data.frame() %>%
      .[-1,] %>%
      mutate(PC = i,
             variable = rownames(.),
             Samples = "NorthAmerica",
             r_sq = tmp_3$r.squared) -> tmp_na
    
    
    Inner_loop = rbind(tmp, tmp_eu, tmp_na)
      
    Inner_loop %<>% 
      mutate(P_bon_test = p.adjust(`Pr(>|t|)`, method = "bonferroni"))
    
    
    lm_inner_list[[i]] = Inner_loop
    
    
  } #i

lm_dat_df = do.call(rbind, lm_inner_list)


lm_dat_df %>%
  filter(P_bon_test < 0.01) %>%
  select(PC, P_bon_test, variable, Samples)

write.table(lm_dat_df, file = "DEST_PC_regressions.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
