### Reproduce panels A and B figure 1
### 

rm(list = ls())


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
library(gmodels)

load("Fig1.panels.AB.dat.Rdata")

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

PCA_table$city %>% unique -> cities_to_select
corr_list = list()
for(i in 1:length(cities_to_select)){
  print(i)
  
  time <- PCA_table$year[which(PCA_table$city == cities_to_select[i])]
  
  
  pc_proj1 <- PCA_table$Dim.1[which(PCA_table$city == cities_to_select[i])]
  tmp1 = cor.test(pc_proj1, time)
  tmp_df1 = data.frame(p_val = tmp1$p.value, 
                       cor = tmp1$estimate,
                       uci = tmp1$conf.int[2],
                       lci = tmp1$conf.int[1],
                       PC = 1,
                       city = cities_to_select[i],
                       sig_sym = ifelse(tmp1$p.value < 0.05, "S", "NS") )
  
  pc_proj2 <- PCA_table$Dim.2[which(PCA_table$city == cities_to_select[i])]
  tmp2 = cor.test(pc_proj2, time)
  tmp_df2 = data.frame(p_val = tmp2$p.value, 
                       cor = tmp2$estimate,
                       uci = tmp1$conf.int[2],
                       lci = tmp1$conf.int[1],
                       PC = 2,
                       city = cities_to_select[i],
                       sig_sym = ifelse(tmp2$p.value < 0.05, "S", "NS") )
  
  pc_proj3 <- PCA_table$Dim.3[which(PCA_table$city == cities_to_select[i])]
  tmp3 = cor.test(pc_proj3, time)
  tmp_df3 = data.frame(p_val = tmp3$p.value, 
                       cor = tmp3$estimate,
                       uci = tmp1$conf.int[2],
                       lci = tmp1$conf.int[1],
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

corr_list_df %>%
  group_by(PC) %>%
  summarise(mean.all = ci(R2)[1],
            lci = ci(R2)[2],
            uci = ci(R2)[3]
             )

label.df_PCA.time %>%
  group_by(city) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(mean_y = mean(year)) %>%
  mutate(endpt = ifelse(year > mean_y, "max", "min")) %>%
  melt(id=c("city","year","mean_y", "endpt")) %>% 
  reshape2::dcast(city ~ variable+endpt ) %>%
  left_join(corr_list_df) ->
  label.df_PCA.time.min.max

#Panel A
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
theme_classic() +
  scale_color_manual(values = c("red","blue")) +
  theme(legend.position = "left")->
  time_projection12

#Panel B
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
            box.padding = 0.5, #max.overlaps = Inf,
            alpha = 1, 
            aes(x=Dim.2, y=Dim.3, 
                label = city))  + 
  geom_path(data = label.df_PCA.time, 
            aes(x=Dim.2, 
                y=Dim.3,
                group = city), arrow = arrow(length=unit(0.10,"cm"),  type = "closed"),size = 0.7) +
  theme_classic() +
  scale_color_manual(values = c("red","blue")) +
  theme(legend.position = "left")->
  time_projection23

#Make plot
ggsave(time_projection12+time_projection23, 
       file = "Fib1.AbadB.pdf",
       width = 9,
       height = 4)

#### Quantification of the signal
#### 

lm(Dim.1 ~ lat + long, weights = MeanEC, data =  PCA_table) -> dim1.lat.long
anova(dim1.lat.long)

lm(Dim.2 ~ lat + long, weights = MeanEC, data =  filter(PCA_table, continent == "Europe")) -> dim2.lat.long
anova(dim2.lat.long)

corr_list_df %>%
  group_by(PC) %>%
  summarize(N = n(),
            mean.cor2 = mean(cor^2))

