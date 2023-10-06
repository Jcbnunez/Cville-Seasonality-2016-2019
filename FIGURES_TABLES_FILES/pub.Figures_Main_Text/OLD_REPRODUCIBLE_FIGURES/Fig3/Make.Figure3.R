#### Reproduce Figure 3
####  
####  
####  
####  Part 1: Reproduce panel A

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
#library(ggVennDiagram)
library(patchwork)
library(data.table)
library(reshape2)
library(doParallel)
library(SNPRelate)
library(SeqArray)
library(vroom)
library(patchwork)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(pegas)
library(vcfR)
library(car)
library(DescTools)
library(ape)
library(ggtree)
library(aplot)
library(forcats)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)

#####
#####
#####
#####

load("./dat.for.panel3A.Rdata")

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -110, ymax = 110), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=dat.for.plot,
    aes(
      x=pos_mean/1e6,
      y=uci.mod,
      color = paste(perm_type, metric)
    ), alpha = 0.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-115,110) +
  scale_color_manual(values = c("grey50", "grey50", "blue", "purple")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Genomic Position (Mb)") +
  ylab("Transformed P-value") +
  xlim(0,20.5) +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) ->
  panel.for.fig3
panel.for.fig3
ggsave(panel.for.fig3, file = "panel.for.fig3.pdf", w = 7, h = 3)

##### 
##### 
##### 
##### 
##### Panel B
##### 
##### 
##### 


#save(sub_pi_d_parsed.plot, inv.dt, final.windows.pos, file = "dat.for.3b.Rdata")
load("./dat.for.3b.Rdata")

ggplot() + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="solid") +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -0, ymax = 0.01), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=sub_pi_d_parsed.plot,
    aes(
      x=BIN_START/1e6,
      y=value,
      color = type),
    alpha = 0.9) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,20.5) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") ->
  pi_d_plot_all

ggsave(pi_d_plot_all, file = "pi_d_plot_all.pdf", w = 7, h = 3.5)

#### Panel C

### Panel 3C
### 
load("FST.data.fig.Rdata")

fst.dat.for.Fig3 %>%
  mutate(mid = (BIN_START+BIN_END)/2) %>%
  melt(id= c("CHROM", "mid", "BIN_START", "BIN_END", "N_VARIANTS", "samp",  "WS")) %>%
  filter(variable == "WEIGHTED_FST") %>%
  ggplot(aes(
    x=mid,
    y=value,
    color = samp,
    linetype = variable
  )) +
  #geom_vline(xintercept = 5192177 , linetype = "solid", color = "blue") +
  geom_vline(xintercept = 2225744, linetype = "dashed") +
  geom_vline(xintercept = 13154180, linetype = "dashed") +
  geom_line(size = 1.0, alpha = 0.9) +
  #facet_grid(WS~samp) +
  xlim(0, 21e6) +
  ggtitle("Weighted FST") +
  theme_bw()

#### Figure D

###
load("./haplowins_pt1.Rdata")
###
rownames(windws_snp_matrix_clean) = windws_snp_matrix_clean$samp_id

analyses_types = list(
  all=c("left",
        "win_3.1",
        "win_4.7",
        "win_5.1",
        "win_6.1",
        "win_6.6",
        "win_9.6",
        "right") #,
  #win3.1=c("win_3.1"),
  #win5.1=c("win_5.1"),
  #win9.6=c("win_9.6")
)

foreach(i=1:length(analyses_types))%do%{
  
  message(names(analyses_types)[i])
  
  analysis = names(analyses_types)[i]
  target_snps <- filter(sim_polarity, win %in% as.character(analyses_types[[i]]) )$SNP_id 
  
  data_in <- windws_snp_matrix_clean[, which(colnames(windws_snp_matrix_clean) %in% target_snps)  ] 
  
  actual_sel_snps <- colnames(data_in)
  
  # dplyr::select(windws_snp_matrix_clean, !(samp_id), filter(SNP_guide_metadata, window == "win5" )$SNP_id )
  D_all <- dist( data_in )
  tre_all <- njs(D_all)
  
  tree_all_plot <- ggtree(tre_all, ignore.negative.edge=TRUE) + 
    #geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    theme_tree2()
  
  tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_missingdat_clean + geom_tippoint(aes(color=pop))
  
  
  is_tip <- tre_all$edge[,2] <= length(tre_all$tip.label)
  ordered_tips <- tre_all$edge[is_tip, 2]
  tre_all$tip.label[ordered_tips]  -> tree_order
  
  joint_figure_polarized_missingdat_clean %>%
    #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
    mutate(samp_id_fct = factor(samp_id, levels = tree_order)) %>% 
    .[complete.cases(.$samp_id_fct),] %>%
    ggplot(
      aes(
        x=as.factor(loci_id),
        y=samp_id_fct,
        fill = polarity
      )
    ) + geom_tile(size = 0.1) +
    facet_grid(.~win, scales = "free", space = "free"
               #ncol = 1, shrink = F
    ) +
    ggtitle(analysis) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())  ->
    karyo_plot_joint_all
  
  karyo_plot_joint_all %>% 
    insert_left(tree_all_plot, width = 0.1)  -> hap_tree_plots_all
  
  ggsave(hap_tree_plots_all, file = paste(analysis, "hap_tree_plots.all.pdf", sep = "."), h = 6, w = 10)
  
}


### Panel E


####

load("DatFor.Haplotypes.trajectory.time.weather.Rdata")

###


#####
Cville_haplotags_for_viz %>%
  ggplot(aes(
    x=yday,
    y=Mean_haplotag,
    #ymin=ci_l,
    #ymax=ci_h,
    color=(temp.max),
  )) + 
  #geom_smooth(method = "lm", se = F, size = 0.8, color = "grey") +
  #geom_errorbar(width = 0.1) +
  scale_color_gradient2(low="steelblue", high = "firebrick2", mid = "gold1", 
                        midpoint = 25) +
  geom_point(aes(shape=as.factor(year))) +
  geom_smooth( se = F, size = 0.8, color = "black", linetype = "dashed") +
  ylim(0,0.38) +
  theme_bw() + 
  facet_grid(win~.)->
  haplo.time.colortemp.ave
haplo.time.colortemp.ave
ggsave(haplo.time.colortemp.ave, file ="haplo.time.colortemp.ave.pdf", h = 6, w = 3.2)


###
###


Cville_haplotags_for_viz %>%
  melt(id = c("sampleId", 
              "collectionDate", 
              "set", 
              "year", 
              "win", 
              "yday", 
              "Mean_haplotag")) %>% 
  filter(variable == "temp.max") %>%
  ggplot(aes(
    x=value,
    y=Mean_haplotag
  )) +
  geom_point(color = "grey",aes(shape=as.factor(year))) +
  geom_smooth(method = "lm", 
              color = "black") +
  theme_bw() +
  facet_grid(win~., scales = "free_x") ->
  eco.vars.afs

ggsave(eco.vars.afs, file ="eco.vars.afs.pdf", h = 6, w = 3.2)


### P-values for regressions:
### 

Cville_haplotags_for_viz %>%
  melt(id = c("sampleId", 
              "collectionDate", 
              "set", 
              "year", 
              "win", 
              "yday", 
              "Mean_haplotag")) %>% 
  filter(variable == "temp.max") %>% 
  nest(data = -c(win) ) %>% 
  mutate(model = map(data, ~lm(Mean_haplotag~value, data = .)),
         tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "value")  %>%
  dplyr::select(win, estimate, p.value, std.error) %>% 
  mutate(estimate = as.numeric(format(estimate, scientific = T, digits = 1)),
         p.value = as.numeric(format(p.value, scientific = T, digits = 1))) %>% 
  as.data.frame()

