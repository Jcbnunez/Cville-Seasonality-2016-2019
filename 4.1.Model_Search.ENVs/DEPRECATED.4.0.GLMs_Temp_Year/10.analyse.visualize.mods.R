## analysis of models

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)


### load glm.out
fl <- list.files("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis", full.names=T)

all.mod.out = foreach(i=1:4, .combine = "rbind")%do%{
  message(fl[i])
  load(fl[i])
  
  win.out %>%
    #mutate(model_name = paste(label,start,end, sep ="_" )) %>%
    group_by(variable, mod,chr,invName,  perm_type, pos_mean , pos_min,  pos_max ) %>%
    summarize(med_rnvp = median(rnp.binom.p),
              rnvp.uci = quantile(rnp.binom.p,0.98),
              rnvp.lci = quantile(rnp.binom.p,0.02),
    ) -> tmp.sum
  tmp.sum
}

### 
save(all.mod.out, file = "all.mod.out.Rdata")
load("all.mod.out.Rdata")
### plot the distributions

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

all.mod.out %>%
  left_join(sets) %>%
  as.data.table() %>%
  mutate(model_name = paste(variable,start,end, sep ="_" )) ->
  dat.for.plot

dat.for.plot %>%
  ggplot(aes(
    x=pos_mean,
    #y=-log10(med_rnvp),
    #ymin=-log10(rnvp10),
    #ymax=-log10(rnvp90),
    #color=perm_type 
  )) +
  #geom_ribbon(width = 0.1) +
  geom_point(data = filter(dat.for.plot, perm_type == "real"),
    aes(y=-log10(rnvp90))) +
  geom_line(data = filter(dat.for.plot, perm_type == "permuted"),
    aes(y=-log10(rnvp90))) +
  facet_grid(model_name~chr, scales = "free_x") ->
  model_properties

ggsave(model_properties, file = "model_properties.pdf")


### version 2
### 
all.mod.out %>%
  filter(chr == "2L") %>%
  filter(label %in% c("F","D","E") ) %>%
  mutate(model_name = paste(label,start,end, sep ="_" )) %>%
  group_by(model_name,chr,invName, perm, perm_type ) %>%
  summarize(med_rnvp = qlogis(median(rnp.binom.p)),
            rnvp90 = qlogis(quantile(rnp.binom.p,0.9)),
            rnvp10 = qlogis(quantile(rnp.binom.p,0.1)),
  ) %>%
  ggplot(aes(
    x=perm,
    #y=qlogis(rnp.binom.p),
    y=med_rnvp,
    ymin=rnvp10,
    ymax=rnvp90,
    color=invName,
    shape = perm_type
  )) +
  geom_errorbar(width = 0.1, position=position_dodge(width=0.5) ) +
  geom_point(position=position_dodge(width=0.5)) +
  #geom_violin(outlier.shape = NA) +
  facet_grid(model_name~invName) ->
  model_properties_2L_perm

ggsave(model_properties_2L_perm, file = "model_properties_2L_perm.pdf")

##### MANHATTAN PLOT
##
all.mod.out %>%
  group_by(perm_type, label, chr, pos_mean) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.01),
            value.wZa = quantile(wZa.p, 0.01)) %>%
  filter(#chr == "2L",
        #perm == 0, 
        label %in% c("F","D","E") ) %>%
  ggplot(aes(
    x=pos_mean/1e6,
    color = perm_type
  ) ) +
  geom_line(aes(y=-log10(value.rnpv)), alpha = 0.6) +
  geom_line(aes(y=log10(value.wZa)), alpha = 0.6) +
  #geom_vline(xintercept = 2225744) +
  #geom_vline(xintercept = 13154180) +
  facet_grid(label~chr, scales = "free_x") ->
  manhatplot.models

ggsave(manhatplot.models, file = "manhatplot.models.pdf", h= 6)

####################
####################
####################
####################
### Outlier Analysis

all.mod.out %>%
  filter(chr == "2L") %>%
  group_by(perm_type, label, chr, pos_mean, pos_min, pos_max) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.0001)) %>%
  filter( label %in% c("F","D","E") ) %>%
  dcast(label+chr+pos_mean+pos_min+pos_max~perm_type) %>% 
  mutate(outlier_test = real <  permuted) %>%
  filter(outlier_test == T) %>%
  group_by(label) %>%
  filter(real < 1e-35) %>%
  mutate(y.art = 1) ->
  outliers_windows.rnpv

all.mod.out %>%
  filter(chr == "2L") %>%
  group_by(perm_type, label, chr, pos_mean, pos_min, pos_max) %>%
  summarize(value.rnpv = quantile(wZa.p, 0.0001)) %>%
  filter( label %in% c("F","D","E") ) %>%
  dcast(label+chr+pos_mean+pos_min+pos_max~perm_type) %>% 
  mutate(outlier_test = real <  permuted) %>%
  filter(outlier_test == T) %>%
  group_by(label) %>%
  filter(real < 1e-35) %>%
  mutate(y.art = 2) ->
  outliers_windows.wZa.p


rbind(data.frame(outliers_windows.rnpv, test="rnpv"), 
data.frame(outliers_windows.wZa.p, test="wza")) ->
  joint.window.analysis


### extract common windows
joint.window.analysis %>%
  dcast(chr+pos_mean+pos_min+pos_max~outlier_test+test+label, 
        value.var = "outlier_test") ->
  summarized.windows.true

summarized.windows.true %>%
  melt(id=c("chr",   "pos_mean",  "pos_min",  "pos_max")) %>% 
  ggplot(aes(
    x=as.factor(pos_mean/1e+6),
    y=variable,
    fill=value
  )) +
    geom_tile() +
    ggtitle("Identifying Peaks across models") +
  theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
        ) ->
    tile.plot

  ggsave(tile.plot, file ="tile.plot.pdf", w = 6, h= 3)

### final windows  
  summarized.windows.true %>%
    .[complete.cases(.),] %>%
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
  
  ggsave(windows_Seg, file = "windows_Seg.pdf", w = 6, h= 3)




