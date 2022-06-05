### Figure Explore GLM outputs
### 
### 

### load libraries
library(patchwork)
library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(data.table)
library(reshape2)
library(car)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)

##### Figure 1
### P.value distributions

output_results_window <- "/scratch/yey2sn/Overwintering_ms/4.GML_plots/final.window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)

#### load windows
#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
load(output_results_window)
###
win.minp.ag <- win.out.final[pr==0.05 & nSNPs>100 & perm!=0 & chr.x == "2L",
                             list(lci=quantile(rnp.binom.p, 0.025, na.rm=T), uci=quantile(rnp.binom.p, .975, na.rm=T), .N),
                             list(mod, chr.x=chr.x, locality, win.i, start, end)] %>%
  filter(locality == "VA_ch", mod=="aveTemp+year_factor")

win.out <- win.out.final %>% filter(locality == "VA_ch", mod=="aveTemp+year_factor", pr == 0.05)
###win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

###win.minp.ag.ag <- win.minp.ag[perm!=0, 
## list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]

### basic MH plot
### Binomial p-value. 
mh.plot.binom <-
  ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="dashed") +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  #geom_line(data=win.out[order(rnp.pr)][nSNPs>100][perm==0],
  #          aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rnp.binom.p)), size=0.2) +
  geom_ribbon(data=win.minp.ag,
              aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.45) +
  geom_point(data=win.out[order(rnp.pr)][nSNPs>100][perm==0],
             aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rnp.binom.p),
                 color=(rnp.pr)), size=1.1) +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab(expression(paste(-log[10], "(Window P)"))) 

ggsave(mh.plot.binom, file="mh.plot.binom.pdf", h=2.9, w=8)


#### DO GLMS beat PERMUTATIONS

win.out.final %>%
  filter(pr == 0.05,
         locality == "VA_ch",
         mod=="aveTemp+year_factor",
         nSNPs > 100)  -> win.out.final.real.perm

win.out.final.real.perm %>%
  ggplot(aes(x = (start.bp+stop.bp)/2,
             y=-1*log10(rnp.binom.p),
             color = win.out.final.real.perm$perm == 0 )) +
  geom_point() +
  facet_wrap(~chr.x, scales = "free_x") ->
  plot_perm_stat

ggsave(plot_perm_stat, file = "plot_perm_stat.pdf")

###  
win.out.final.real.perm %>%
  filter(perm != 0) %>%
  group_by(chr.x, start.bp, stop.bp ) %>%
  summarize(quant_h=quantile(rnp.binom.p, 0.001),
            quant_l=quantile(rnp.binom.p, 0.999)
              ) ->
  win.out_quant.summ.perm

left_join( filter(win.out.final.real.perm, perm == 0 ), 
           win.out_quant.summ.perm,
           by = c("chr.x","start.bp", "stop.bp")) %>% 
  mutate(pos_med = (start.bp+stop.bp)/2 ) ->
  win.out_quant.summ.processed
  
##win.out_quant.summ.processed %>%
##  dcast(chr.x+start.bp+stop.bp~perm_stat, value.var = "quant_p") %>%
##  mutate(pos_med = (start.bp+stop.bp)/2 ) ->
##  win.out_quant.summ.dcast
findpeaks(-log10(filter(win.out_quant.summ.processed, rnp.binom.p <  quant_h)$rnp.binom.p) ) -> peaks.glm
peaks.glm %<>%
  as.data.frame() %>%
  dplyr::select(peak_val = V1, win.out.id = V2, peak_begin = V3, peak_end = V4)

filter(win.out_quant.summ.processed, rnp.binom.p <  quant_h) %>%
  mutate(win.out.id=1:dim(.)[1]) %>%
  left_join(peaks.glm) %>%
  mutate(is.peak = case_when(is.na(peak_val) ~ "F",
                             !is.na(peak_val) ~ "T"
                             ))-> beat_perm_wins

ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=stop), linetype="dashed") +
  geom_line(data = win.out_quant.summ.processed, 
            aes(x=pos_med, y= -log10(rnp.binom.p)), 
            color = "red", alpha = 0.9) +
  geom_ribbon(data = win.out_quant.summ.processed, 
              aes(x=pos_med, ymax= -log10(quant_h), ymin = -log10(quant_l),), 
              fill = "purple" , alpha = 0.4) +
  geom_point(data = win.out_quant.summ.processed, 
             aes(x=pos_med, y= -log10(rnp.binom.p)), 
             color = "orange", alpha = 0.3, size = 1.1) +
  geom_point(data = beat_perm_wins,  
             aes(x=pos_med, y= -log10(rnp.binom.p), fill = is.peak), 
              alpha = 0.9, shape = 23, size = 1.1) +
  theme_bw() +
  facet_wrap(~chr.x, scales = "free_x") ->
  plot_perm_real_final

ggsave(plot_perm_real_final, file = "plot_perm_real_final.pdf", h = 4, w = 8)



#### WzA
#### 
mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="dashed") +
  #geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_line(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rnp.wZa.p)), size=0.8, color = "steelblue") +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab(expression(paste(-log[10], "(wZa P)"))) 

ggsave(mh.plot.wza, file="mh.plot.wza.pdf", h=2.9, w=5)


#### Join analysis

#### find peaks for the wZa test
#### 
wZa_dat = -log10(win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0]$rnp.wZa.p)
binom_dat = -log10(win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0]$rnp.binom.p)


findpeaks(wZa_dat, #minpeakheight = mean(binom_dat)+5*sd(binom_dat), 
          minpeakdistance = 3) -> peaks.wza
findpeaks(binom_dat, #minpeakheight = mean(binom_dat)+5*sd(binom_dat), 
          minpeakdistance = 3 ) -> peaks.binom

peaks.wza %<>%
  as.data.frame() %>%
  dplyr::select(win.out.id = V2, 
                ) %>%
  mutate( status.wza = "peak") 

peaks.binom %<>%
  as.data.frame() %>%
  dplyr::select( 
    win.out.id = V2, 
  ) %>%
  mutate( status.binom = "peak")

win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0] %>%
  mutate(win.out.id=1:dim(.)[1]) %>%
  left_join(peaks.binom) %>%
  left_join(peaks.wza) %>% 
  mutate(joint_peak = case_when(status.binom == "peak" & status.wza == "peak" ~ TRUE,
                                is.na(status.binom) | is.na(status.wza)  ~ FALSE)) %>%
  left_join(  win.out_quant.summ.perm,   #object generated above
             by = c("chr.x","start.bp", "stop.bp")) %>% 
  mutate(pos_med = (start.bp+stop.bp)/2 ) ->
  joint.analysis.perm.real.processed

save(joint.analysis.perm.real.processed, inv.dt,  file = "joint.analysis.perm.real.processed.Rdata")
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
#### Plot the joint figure
load("./joint.analysis.perm.real.processed.Rdata")


joint.analysis.perm.real.processed.2L = filter(joint.analysis.perm.real.processed, chr.x == "2L")
  
ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=stop), linetype="dashed") +
  geom_line(data = joint.analysis.perm.real.processed.2L, 
            aes(x=pos_med, y= -log10(rnp.binom.p)), 
            color = "red", alpha = 0.9) +
  geom_ribbon(data = joint.analysis.perm.real.processed.2L, 
            aes(x=pos_med, ymin= log10(rnp.wZa.p)), ymax = 0,
            fill = "blue", alpha = 0.4) +
  geom_ribbon(data = joint.analysis.perm.real.processed.2L, 
              aes(x=pos_med, ymax= -log10(quant_h), ymin = -log10(quant_l),), 
              fill = "purple" , alpha = 0.4) +
  #geom_point(data = joint.analysis.perm.real.processed.2L, 
  #           aes(x=pos_med, y= -log10(rnp.binom.p)), 
  #           color = "orange", alpha = 0.3, size = 1.1) +
  geom_point(data = filter(joint.analysis.perm.real.processed.2L, rnp.binom.p < quant_h ),  
             aes(x=pos_med, y= -log10(rnp.binom.p)), 
             alpha = 0.9, shape = 23, size = 0.7) +
  theme_bw() +
  facet_wrap(~chr.x, scales = "free_x") ->
  plot_perm_real_final_with_wZA

ggsave(plot_perm_real_final_with_wZA, file = "plot_perm_real_final_with_wZA.pdf", h = 4, w = 8)




