### Make Figure 1

rm(list = ls())

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
library(scales)
library(gmodels)
library(forcats)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(gmodels)
library(MASS)


####
load("Fig1.panels.AB.dat.Rdata")

label.df_PCA <- PCA_table %>% 
  group_by(city) %>% 
  dplyr::summarize(Dim.1 = mean(Dim.1), 
            Dim.2 = mean(Dim.2),
            Dim.3 = mean(Dim.3)) 

label.df_PCA.time <- PCA_table %>% 
  group_by(city, year) %>% 
  dplyr::summarize(Dim.1 = mean(Dim.1), 
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
  reshape2::melt(id=c("city","year","mean_y", "endpt")) %>% 
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


#####
#####

### panel DC
load("Fig1.panelC.Rdata")

Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  group_by(pop1, bin_date) %>% 
  summarise(FST_mean = mean(FST)) %>%
  dcast(pop1 ~ bin_date) ->
  mean_fst

Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  left_join(mean_fst) %>% 
  ggplot(
    aes(
      x=fct_reorder(pop1, `1.within`),
      y=FST,
      color=bin_date
    )
  ) +
  geom_boxplot(width = 0.7) +
  xlab("Sampling Locale") +
  coord_flip() +
  theme_bw() +
  ylab(expression(italic(F)[ST])) -> 			
  fst_boxplot
fst_boxplot
ggsave(fst_boxplot,
       file = "fst_boxplot.pdf",
       width = 5,
       height = 4)

### Quantify

Out_comp_vector_samepops %>%
  group_by(bin_date) %>%
  summarise(med = quantile(FST, 0.5 ))

anova(lm(FST ~ bin_date, data = filter(Out_comp_vector_samepops, bin_date %in% c("1.within", "2.Overwinter") )))

### panel D
load("Fig1.panelD.Rdata")

multiyear_samps %>%  
  group_by(pop1, year_diff) %>%
  dplyr::summarize(Mean = mean(FST),
                   SD = sd(FST)) %>%
  ggplot(
    aes(
      x=(year_diff),
      y=(Mean),
      ymin=(Mean)-SD,
      ymax=(Mean)+SD,
      fill=pop1,
      color=pop1,
      #color=as.factor(day_diff_year_scaled)
    ))  + 
  #geom_boxplot(width = 0.4) +
  geom_errorbar(width = 0.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_point(shape = 21, size = 2.3, position=position_dodge(width=0.5)) +
  xlab("Number of Winters") +
  theme_bw() +
  xlab(expression(Delta[Years])) +
  ylab(expression(F[ST])) +
  #scale_shape_manual(values = c(21,22, 23)) +
  scale_color_brewer(palette ="Dark2") +
  scale_fill_brewer(palette ="Dark2") ->
  fst_allpop_overwint
fst_allpop_overwint

ggsave(fst_allpop_overwint,
       file = "fst_allpop_overwint.pdf",
       width = 6,
       height = 2.3)

## Quantify:
anova(lm(FST ~ year_diff*pop1, data = multiyear_samps ))

