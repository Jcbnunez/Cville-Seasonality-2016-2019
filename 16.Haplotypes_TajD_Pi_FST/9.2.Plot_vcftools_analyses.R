### Collect VCFTools analyses
### 
rm(list = ls())


library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)

inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

####### Load all
final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/0.OLD_ANALYSIS.preAprl1/pi_D_datforplot.Rdata")

sub_pi_d_parsed %>%
  filter(pop == "cm" & resolution == "W_100000" & type != "all") ->
  sub_pi_d_parsed.plot
  
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


save(sub_pi_d_parsed.plot, inv.dt, final.windows.pos, file = "dat.for.3b.Rdata")







sub_pi_d_parsed %>% 
  group_by(pop, type, resolution) %>%
  summarize(N = n())

sub_pi_d_parsed %>%
  filter(resolution == "W_100000") %>%
  ggplot(
    aes(
      x=BIN_START,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="solid") +
  geom_line(alpha = 0.5) +
  ggtitle("Samples split by Karyotype") +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  pi_d_plot_sub

ggsave(pi_d_plot_sub, file = "pi_d_plot_sub.pdf")

###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
#####  ######
load("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/pi_D_datforplot.Rdata")

sub_pi_d_parsed %>%
  filter(resolution == "W_10000") %>% 
  mutate(zoom_win = case_when(BIN_START > 4.8e6 & BIN_START < 5.7e6 ~ "1.win5",
                              BIN_START > 5.8e6 & BIN_START < 7e6 ~ "2.win6",
                              BIN_START > 9.0e6 & BIN_START < 10.6e6 ~ "3.win10",
                              )) %>% 
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) -> D_data_zoom

D_data_zoom %>%
  ggplot(
      aes(
      x=BIN_START/1e6,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_line(alpha = 0.5) +
  ggtitle("Zoom-in to GLM peaks") +
  theme_bw() +
  facet_wrap(variable~zoom_win, ncol = 3, scales = "free") ->
  pi_d_plot_zoom

ggsave(pi_d_plot_zoom, file = "pi_d_plot_zoom.pdf", h= 6, w = 9)

#### ADD GLM P_VALUES
output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")
ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)
#### load windows
#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
load(output_results_window)
###
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0 & chr.x == "2L",
                       list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
                       list(mod, chr.x=chr.x, locality, win.i, start, end)] %>%
  filter(locality == "VA_ch",mod=="aveTemp+year_factor") %>%
  filter(mod == "aveTemp+year_factor") %>%
  mutate(zoom_win = case_when(
    start > 4.8e6 & end < 5.7e6 ~ "1.win5",
    start > 5.8e6 & end < 7e6 ~ "2.win6",
    start > 9.0e6 & end < 10.6e6 ~ "3.win10")) %>% 
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) %>%
  mutate( variable = "GLM" ) 

win.minp.ag.zoom = win.minp.ag %>% filter(zoom_win %in% c("1.win5", "2.win6","3.win10" ))

win.out.zoom =  win.out %>% filter(locality == "VA_ch", mod=="aveTemp+year_factor", chr.x == "2L") %>%
  mutate(zoom_win = case_when(
    start > 4.8e6 & end < 5.7e6 ~ "1.win5",
    start > 5.8e6 & end < 7e6 ~ "2.win6",
    start > 9.0e6 & end < 10.6e6 ~ "3.win10")) %>% 
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) %>%
  mutate( variable = "GLM" ) 

win.out.zoom %>% 
  group_by(zoom_win) %>%
  slice_max(order_by = -1*log10(rbinom.p) ) -> 
  max_rnpvls  ## <---- classical approach

classic_way_max_rnpvls = max_rnpvls

ggplot() +
  geom_line(data=filter(D_data_zoom, type == "inv"),
    aes(
      x=BIN_START/1e6,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_vline(data = dplyr::select(max_rnpvls, -variable) , aes(xintercept = start/1e6), linetype = "dashed") +
  geom_ribbon(data=win.minp.ag.zoom,
              aes(#x=(start/2 +end/2)/1e6, 
                  x=start/1e6,
                  ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.45) +
  geom_point(data=win.out.zoom[order(rnp.pr)][nSNPs>100][perm==0],
             aes(#x=(start/2 +end/2)/1e6 
                  x=start/1e6,
                  y=-1*log10(rbinom.p),
                 #color=(rnp.pr)
                 ), size=1.1) +  
  ggtitle("Zoom-in to GLM peaks") +
  theme_bw() +
  facet_wrap(variable~zoom_win, ncol = 3, scales = "free") ->
  GLM_d_plot_zoom

ggsave(GLM_d_plot_zoom, file = "GLM_d_plot_zoom.pdf", h= 6, w = 9)

####### Alan's special case
###
###
###
###
###
alans_wza_win <- get(load("./WZA_window.dt.VA_ch_0.Rdata"))
alans_wza_win %>%
  filter(mod == "aveTemp+year_factor",
         chr.x == "2L",
         pr == 0.05) %>% mutate(zoom_win = case_when(
  start > 4.8e6 & end < 5.7e6 ~ "1.win5",
  start > 5.8e6 & end < 7e6 ~ "2.win6",
  start > 9.0e6 & end < 10.6e6 ~ "3.win10")) %>% 
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) ->
  alans_wza_win_parsed

alans_wza_win_parsed %>% 
  group_by(zoom_win) %>%
  slice_max(order_by = -1*log10(rbinom.p) ) -> 
  max_rnpvls ## <---- Alan's mew way approach

AlansNew_way_max_rnpvls = max_rnpvls

ggplot() +
  geom_line(data=mutate(filter(D_data_zoom, type == "all"), variable = paste("1.all", variable, sep = "_")),
            aes(
              x=BIN_START/1e6,
              y=value,
              color =type,
              linetype = pop)
  ) + 
  geom_line(data=mutate(filter(D_data_zoom, type == "inv"), variable = paste("2.inv", variable, sep = "_")),
            aes(
              x=BIN_START/1e6,
              y=value,
              color =type,
              linetype = pop)
  ) + 
  geom_line(data=mutate(filter(D_data_zoom, type == "std"), variable = paste("3.std", variable, sep = "_")),
            aes(
              x=BIN_START/1e6,
              y=value,
              color =type,
              linetype = pop)
  ) + 
  geom_point(data=mutate(alans_wza_win_parsed, variable = "GLM"),
             aes(#x=(start/2 +end/2)/1e6 
               x=start/1e6,
               y=-1*log10(rbinom.p),
             ), size=1.1, alpha = 0.5) +  
  geom_point(data=mutate(alans_wza_win_parsed, variable = "wZa_p"),
             aes(#x=(start/2 +end/2)/1e6 
               x=start/1e6,
               y=-1*log10(wZa.p),
             ), size=1.1) +  
  geom_vline(data = max_rnpvls, aes(xintercept = start/1e6), linetype = "dashed", color = "blue") +
  ggtitle("Zoom-in to GLM peaks | Alan's new set" ) +
  theme_bw() +
  facet_wrap(variable~zoom_win, ncol = 3, scales = "free") ->
  GLM_d_plot_zoom

ggsave(GLM_d_plot_zoom, file = "GLM_d_plot_zoom.pdf", h= 12, w = 9)

### ### ###

###
##windows_of_interest =
##data.frame(rbind(
##  mutate(AlansNew_way_max_rnpvls, method = "fixed_snp_win", variable = "GLM"),
##  mutate(classic_way_max_rnpvls, method = "fixed_distance_win")))
##save(windows_of_interest, file = "windows_of_interest.Rdata")
##
##windows_of_interest %>%
##  ggplot(aes(
##    x=start,
##    y=1:6,
##    xend=end,
##    yend=1:6,
##    color = method
##  )) +
##  geom_segment() +
##  facet_wrap(~zoom_win, ncol = 1, scales = "free" ) ->
##  win_placement
##
##ggsave(win_placement, file ="win_placement.pdf")
##
##

###
###
###
###
###
### add genes regions 
#bring in the snp data:
glm_ld_outliers <- "/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt"
outlier_list <- fread(glm_ld_outliers)

outlier_list %>%
  filter(Transcript_biotype == "protein_coding") %>%
  group_by(Gene_Name) %>%
  summarise(min_pos = min(pos),
            max_pos = max(pos))  %>%
  mutate(zoom_win = 
           case_when(min_pos > 4.7e6 & max_pos < 5.3e6 ~ "1.win5",
                     min_pos > 5.8e6 & max_pos < 7e6 ~ "2.win6",
                     min_pos > 9.0e6 & max_pos < 10.6e6 ~ "3.win10")) -> gene_list

gene_list  %>% 
  .[-grep("^CG" ,.$Gene_Name),] %>%
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) %>%
  .[order(.$min_pos),] %>%
  mutate(gene_id = 1:dim(.)[1],
         analysis = "2.GeneId") ->
  gene_list_select

D_data_zoom %<>%
  mutate(analysis = "1.TajD")

data.frame(xint = c(6.07, 6.69, 9.64, 10.0),
           zoom_win = c("2.win6", "2.win6", "3.win10","3.win10")) ->
  vlines

ggplot() +
  geom_vline(data=vlines, aes(xintercept=xint), linetype="dashed") +
  geom_line(data= filter(D_data_zoom, variable == "TajimaD", pop == "cm"),
    aes(
      x=BIN_START/1e6,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_segment(data=gene_list_select,
               aes(
                x=min_pos/1e6,
                xend=max_pos/1e6,
                y=gene_id/1e6,
                yend=gene_id/1e6), alpha = 0.5
               )  +
  geom_text(data=gene_list_select,
               aes(
                 x=(max_pos)/1e6,
                 y=(gene_id+0.5)/1e6,
                 label=Gene_Name
                 ),size = 1
  )  +
  facet_wrap(analysis~zoom_win, ncol = 3, scales = "free") ->
  gene_d_plot_zoom

ggsave(gene_d_plot_zoom, file = "gene_d_plot_zoom.pdf", h= 6, w = 9)

  

