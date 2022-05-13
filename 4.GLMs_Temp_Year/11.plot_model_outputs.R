### Make new figure 3
### 

### libraries
library(tidyverse)
library(data.table)
library(tidyr)
library(viridis)
library(patchwork)
library(magrittr)
library(RColorBrewer)
library(rstatix)

##### input
output_results_window <- "/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/all.mod.out.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

names(inv.dt)[1] = "chr"
inv.dt %>%
  melt(id = c("chr", "invName"), value.name = "pos_mean") ->
  inv.dt.melt

#####
outlier_haplowins = 
  data.frame(win.name = c("w4.7", "w5.1", "w6.2", "w6.8", "w9.5" ),
             start = c(4656622, 5105919, 6155931, 6805798, 9505855 ),
             end = c(4805715, 5255685, 6355509, 6905746, 9605419))

####
####
load(output_results_window)

### summarize and plot
all.mod.out %>%
  filter(chr == "2L",
    #perm == 0, 
    label %in% c("F","D","E") ) %>%
  mutate(model_name = paste(label,start,end, sep ="_" )) %>%
  group_by(perm_type, model_name, label, chr, pos_mean) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.01),
            value.wZa = quantile(wZa.p, 0.01)) ->
  sum_dat_for_plot

  ggplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 35, linetype = "dashed") +
  geom_hline(yintercept = -35, linetype = "dashed") +
  geom_vline(data=inv.dt[chr.x == "2L"], aes(xintercept=start/1e6, #linetype=invName
                                             )) +
  geom_vline(data=inv.dt[chr.x == "2L"], aes(xintercept=stop/1e6, #linetype=invName
                                             )) +
  geom_rect(data = outlier_haplowins, 
            aes(xmin=start/1e6, xmax=end/1e6, ymin= -200, ymax = 100), 
               fill = "gold", alpha = 0.8) +
  geom_ribbon(data = sum_dat_for_plot,
              aes(x=pos_mean/1e6, 
                ymax=-log10(value.rnpv), 
                ymin=0,
                fill = perm_type),
                alpha = 0.6) +
  geom_ribbon(data = sum_dat_for_plot,
                  aes(x=pos_mean/1e6, 
                  ymax=log10(value.wZa), 
                  ymin=0,
                  fill = perm_type),
              alpha = 0.6) +
  scale_fill_manual(values=c("gray39", "red")) +
  facet_grid(model_name~.) +
  theme_bw()->
  model_plot

ggsave(model_plot, file = "model_plot.pdf", h = 5, w = 7)
  
  

### summarize and plot for all genome
all.mod.out %>%
  filter(#chr == "2L",
         #perm == 0, 
         label %in% c("F","D","E") ) %>%
  mutate(model_name = paste(label,start,end, sep ="_" )) %>%
  group_by(perm_type, model_name, label, chr, pos_mean) %>%
  summarize(value.rnpv = quantile(rnp.binom.p, 0.01),
            value.wZa = quantile(wZa.p, 0.01)) ->
  sum_dat_for_plot_allG

ggplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 35, linetype = "dashed") +
  geom_hline(yintercept = -35, linetype = "dashed") +
  #geom_vline(data=inv.dt[chr.x == "2L"], aes(xintercept=start/1e6, #linetype=invName
  #)) +
  geom_vline(data=inv.dt.melt, aes(xintercept=pos_mean/1e6, linetype =invName  #linetype=invName
  )) +
  #geom_rect(data = outlier_haplowins, 
  #          aes(xmin=start/1e6, xmax=end/1e6, ymin= -200, ymax = 100), 
  #          fill = "gold", alpha = 0.8) +
  geom_ribbon(data = sum_dat_for_plot_allG,
              aes(x=pos_mean/1e6, 
                  ymax=-log10(value.rnpv), 
                  ymin=0,
                  fill = perm_type),
              alpha = 0.6) +
  geom_ribbon(data = sum_dat_for_plot_allG,
              aes(x=pos_mean/1e6, 
                  ymax=log10(value.wZa), 
                  ymin=0,
                  fill = perm_type),
              alpha = 0.6) +
  scale_fill_manual(values=c("gray39", "red")) +
  facet_grid(model_name~chr, scales = "free_x") +
  theme_bw()->
  model_plot_allGenome

ggsave(model_plot_allGenome, file = "model_plot_allGenome.pdf", h = 6, w = 8)



  
  
