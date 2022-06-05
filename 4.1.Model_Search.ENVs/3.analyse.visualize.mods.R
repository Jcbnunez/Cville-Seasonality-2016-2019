## analysis of models
rm(list = ls())

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(magrittr)
library(ggVennDiagram)
library(patchwork)

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

chr.sets = list(all=c("2L", "2R", "3L", "3R"),
                ch2L = "2L",
                ch2R = "2R",
                ch3L = "3L",
                ch3R = "3R")

##load model files
mods_fin <- fread("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/final_models.txt")

mods_fin$mod_var -> models


### plot eco variables
load("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/weather.Rdata")

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(V1, temp.max) -> tmax

weather.ave %>%
  filter(mod == 3) %>%
  dplyr::select(V1, precip.var) -> pvar

weather.ave %>%
  filter(mod == 4) %>%
  dplyr::select(V1, temp.ave) -> tave

weather.ave %>%
  filter(mod == 11) %>%
  dplyr::select(V1, precip.ave) -> pave

left_join(tmax, pvar) %>% left_join(tave) %>% left_join(pave) -> obs_for_cors

obs_for_cors %>%
  ggplot(aes(
    x=temp.max,
    y=precip.var,
  )) +
  geom_point() +
  geom_smooth(method = "lm") ->
  t.prec.var

obs_for_cors %>%
  ggplot(aes(
    x=temp.max,
    y=precip.ave,
  ))+
  geom_point() +
  geom_smooth(method = "lm") ->
  t.prec.ave

ggsave(t.prec.var+t.prec.ave, file = "corrs.raw.pts.pdf")
### load and extract real data

##### base files
#base <- "/project/berglandlab/alan/environmental_ombibus"
#
#real.dat.models = foreach(k = 1:4, .combine = "rbind")%do%{
#  
#  file <- paste(base, models[k], paste(models[k],"glmRNP.Rdata", sep = ".") , sep = "/" )
#  print(file)
#  
#  out.glm <- get(load(file))
# 
#   out.glm %>%
#    filter(perm == 0) ->
#    real_data_mod
#  
#  real_data_mod
#}
#
#save(real.dat.models, file = "real.dat.models.Rdata")
load("real.dat.models.Rdata")

real.dat.models %<>%
   left_join(sets) %>%
    mutate(model_name = paste(variable,start,end, sep ="_" )) %>%
    mutate(SNP_id = paste(chr,pos, sep = "_"))
##
##foreach(i=1:5)%do%{
##  
##  list_outlier = list(
##    t.max = filter(real.dat.models, rnp < 0.05, model_name == "temp.max_0_15", chr %in% chr.sets[[i]] )$SNP_id,
##    pcip.var = filter(real.dat.models, rnp < 0.05, model_name == "precip.var_0_30", chr %in% chr.sets[[i]])$SNP_id,
##    t.ave = filter(real.dat.models, rnp < 0.05, model_name == "temp.ave_7_15", chr %in% chr.sets[[i]])$SNP_id,
##    pcip.ave = filter(real.dat.models, rnp < 0.05, model_name == "precip.ave_0_90", chr %in% chr.sets[[i]])$SNP_id
##  )
##  
##  ggVennDiagram(list_outlier, label_alpha = 0) +
##    scale_fill_gradient(low = "white", high = "red") +
##    ggtitle(names(chr.sets)[i]) -> ggvenn.plot
##  ggsave(ggvenn.plot, file = paste(names(chr.sets)[i],"ggvenn.plot.pdf", sep = ".") )
##  
##}
##
#dcast(chr+pos~model_name, value.var = "rnp")

  
#### load glm.out
#fl <- list.files("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis", full.names=T)
#
#all.mod.out = foreach(i=1:4, .combine = "rbind")%do%{
#  message(fl[i])
#  load(fl[i])
#  
#  win.out 
#}

#all.mod.out %<>%
#  left_join(sets) %>%
#  mutate(model_name = paste(variable,start,end, sep ="_" )) 
#
#### 
#save(all.mod.out, file = "all.mod.out.Rdata")
load("all.mod.out.Rdata")

####################
####################
####################
####################
### Outlier Analysis

all.mod.out %>%
  filter(chr == "2L", pos_min > 2225744 &  pos_max < 13154180) %>%
  group_by(perm_type, model_name, chr, pos_mean, pos_min, pos_max) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.02)) %>%
  dcast(model_name+chr+pos_mean+pos_min+pos_max~perm_type) %>% 
  mutate(outlier_test = real <  permuted) %>%
  filter(outlier_test == T) %>%
  group_by(model_name) %>%
  filter(real < quantile(all.mod.out$rnp.binom.p, 0.001) ) %>%
  mutate(y.art = 1) ->
  outliers_windows.rnpv


all.mod.out %>%
  filter(chr == "2L", pos_min > 2225744 &  pos_max < 13154180) %>%
  group_by(perm_type, model_name, chr, pos_mean, pos_min, pos_max) %>%
  summarize(value.wZa = quantile(wZa.p, 0.02)) %>%
  dcast(model_name+chr+pos_mean+pos_min+pos_max~perm_type) %>% 
  mutate(outlier_test = real <  permuted) %>%
  filter(outlier_test == T) %>%
  filter(real < quantile(all.mod.out$wZa.p, 0.001) ) %>%
  mutate(y.art = 2) ->
  outliers_windows.wZa.p

####

rbind(
data.frame(outliers_windows.rnpv, test="rnpv"), 
data.frame(outliers_windows.wZa.p, test="wza")) ->
  joint.window.analysis

### extract common windows
joint.window.analysis %>%
  dcast(chr+pos_mean+pos_min+pos_max~outlier_test+test+model_name, 
        value.var = "outlier_test") ->
  summarized.windows.true

names(summarized.windows.true) %>% .[5:11] -> vars_select
names(summarized.windows.true) %>% .[1:4] -> vars_ids

summarized.windows.true %>%
  melt(id = vars_ids ) %>%
  mutate(value.score = case_when(is.na(value) ~ 0,
                                 TRUE ~ 1)) %>%
  group_by(pos_mean, pos_min, pos_max) %>%
  summarise(Score.add = sum(value.score) ) %>%
  filter(Score.add >= 6) ->
  summarized.windows.true.final


##summarized.windows.true.final %>%
##  melt(id=c("chr",   "pos_mean",  "pos_min",  "pos_max")) %>% 
##  ggplot(aes(
##    x=as.factor(pos_mean/1e+6),
##    y=variable,
##    fill=value
##  )) +
##    geom_tile() +
##    ggtitle("Identifying Peaks across models") +
##  theme(#axis.title.x=element_blank(),
##        axis.text.x=element_blank(),
##        #axis.ticks.x=element_blank()
##        ) ->
##    tile.plot
##
##  ggsave(tile.plot, file ="tile.plot.pdf", w = 6, h= 3)
##
### final windows  
summarized.windows.true.final%>%
    mutate(y.art = 1:n()) %>%
    ggplot(aes(
      x=pos_min,
      xend=pos_max,
      y=1,
      yend=1,
    )) +
    geom_segment() +
    geom_vline(xintercept = 2225744) +
    geom_vline(xintercept = 13154180) +
    ylim(0,2)  ->
    windows_Seg
ggsave(windows_Seg, file = "windows_Seg.pdf" )

### final windows
data.frame(
  start= c(4650065, 5100324, 6100321, 9500286),
  end= c(4799922, 5349218, 6349489, 9700005)
) %>% mutate(mid.point = (start/2) + (end/2) ) %>%
  mutate(win.name = paste("win", round(mid.point/1e6, 1), sep = "."  )) ->
  final.windows.pos

save(final.windows.pos, file = "final.windows.pos.Rdata")

#### PLOTs

load("all.mod.out.Rdata")
load("final.windows.pos.Rdata")

##### MANHATTAN PLOT
##
all.mod.out %>%
  filter(chr == "2L",
         model_name %in% c("temp.ave_7_15", "temp.max_0_15")  ) %>%
  group_by(perm_type, model_name, chr, pos_mean) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.02),
            value.wZa = quantile(wZa.p, 0.02)) -> dat.in.plot

dat.in.plot$perm_type = factor(dat.in.plot$perm_type, levels = c("real", "permuted"))
library(viridis)
library("ggsci")

  ggplot() +
    geom_rect(data = final.windows.pos,
              aes(xmin=start/1e6, xmax = end/1e6,
                  ymin = -230, ymax = 120), 
              alpha = 0.7, fill = "gold") +
    geom_ribbon(data = filter(dat.in.plot, perm_type == "permuted"),
                aes(
                  x=pos_mean/1e6,
                  ymin = 0,
                  ymax=-log10(value.rnpv),
                  #linetype = perm_type,
                  fill = model_name
                ),  alpha = 0.5) +
    geom_ribbon(data = filter(dat.in.plot, perm_type == "permuted"),
                aes(
                  x=pos_mean/1e6,
                  ymax = 0,
                  ymin=log10(value.wZa),
                  #linetype = perm_type,
                  fill = model_name
                ),  alpha = 0.5) +
  geom_line(data = filter(dat.in.plot, perm_type == "real"),
            aes(
              x=pos_mean/1e6,
              y=-log10(value.rnpv),
              linetype = model_name,
              color = model_name
            ),  alpha = 0.9) +
  geom_line(data = filter(dat.in.plot, perm_type == "real"),
            aes(
              x=pos_mean/1e6,
              y=log10(value.wZa),
              linetype = model_name,
              color = model_name
            ) , alpha = 0.9) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Greys") +
    facet_grid(.~chr, scales = "free_x") ->
  manhatplot.models

ggsave(manhatplot.models, file = "manhatplot.models.pdf", w= 6, h = 3)

#### LD PLOTS
load("../7.LD/ld_df_winanot_pair.Rdata")
PEAKS_for_ANALYSIS = final.windows.pos

ggplot() +
  geom_point(data=filter(ld_df_winanot_pair, R2 >= 0.68 & BP_A < BP_B),
             aes(
               x=BP_A/1e6, y=BP_B/1e6,
               #shape = as.factor(Win_stat_A),
               color=R2, #size = R2
             ),size = 0.09, alpha = 0.1) +
  
  #+ geom_point( shape = 15, alpha = 0.5) +
  #geom_abline(intercept = 0, slope = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "purple", midpoint =  0.80) +
  geom_vline(data = PEAKS_for_ANALYSIS, 
             aes(xintercept=mid.point/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.5, size = 3, color = "gold1"
  ) +
  geom_hline(data = PEAKS_for_ANALYSIS, 
             aes(yintercept=mid.point/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.5, size = 3, color = "gold1"
  ) +
  theme_bw() +
  theme(#legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),) ->
  #facet_wrap(~winComp, scale ="free") ->
  ld_triag_smorg

ggsave(ld_triag_smorg, file = "ld_triag_smorg.pdf", w=4, h=3)

dat.high.ld = filter(ld_df_winanot_pair, R2 >= 0.60 & BP_A < BP_B) 
runif(10000, 1, dim(dat.high.ld)[1]) -> rand.vals
dat.rand = dat.high.ld %>%
              .[rand.vals,]


ggplot() +
  geom_point(data=dat.rand,
             aes(
               x=BP_A/1e6, y=BP_B/1e6,
               #shape = as.factor(Win_stat_A),
               color=R2, #size = R2
             ),size = 0.3, alpha = 0.5) +
  
  #+ geom_point( shape = 15, alpha = 0.5) +
  #geom_abline(intercept = 0, slope = 1) +
  scale_color_gradient2(low = "grey", high = "red", mid = "purple", midpoint =  0.80) +
  geom_rect(data = PEAKS_for_ANALYSIS, 
             aes(xmin=start/1e6,
                 xmax=end/1e6,
                 ymin=start/1e6,
                 ymax=end/1e6,
                 ), 
             #linetype="solid", color = "purple", 
             alpha = 0.5, size = 3, color = "gold1"
  ) +
  ##geom_hline(data = PEAKS_for_ANALYSIS, 
  #           aes(yintercept=mid.point/1e6), 
  #           #linetype="solid", color = "purple", 
  #           alpha = 0.5, size = 3, color = "gold1"
  #) +
  theme_bw() +
  theme(#legend.position = "none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),) ->
  #facet_wrap(~winComp, scale ="free") ->
  ld_triag_smorg.rand

ggsave(ld_triag_smorg.rand, file = "ld_triag_smorg.rand.pdf", w=4, h=3)


