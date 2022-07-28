
### load libraries
library(patchwork)
library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(car)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(DescTools)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)

#### Import all the simulation dta
##### input
output_results_window <- "/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/all.mod.out.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

#####
outlier_haplowins = 
  data.frame(win.name = c("w4.6", "w5.1", "w6.2", "w6.8", "w9.5" ),
             start = c(4656622, 5105919, 6155931, 6805798, 9505855 ),
             end = c(4805715, 5255685, 6355509, 6905746, 9605419))

####
####
load(output_results_window)

all.mod.out %>%
  filter(#chr == "2L",
         #perm == 0, 
         label %in% c(
                      #"F","D",
                      "E") ) %>%
  mutate(in_lab = case_when(invName != "noInv" ~ "Inv",
                            TRUE ~ "noInv")) %>%
  mutate(model_name = paste(label,start,end, sep ="_" )) %>%
  .[,c("perm_type", "perm", "pos_mean", "model_name", "label", "chr",  "in_lab", "wZa.p", "rnp.binom.p")] %>%
  melt(id = c("perm_type", "perm", "pos_mean" , "model_name", "label", "chr",  "in_lab")) %>% 
  group_by(perm_type, perm, model_name, label, chr,  in_lab, variable) %>%
  summarize(median.val = median(value),
            upper.95  = quantile(value, 0.95),
            lower.05  = quantile(value, 0.05)) %>% 
  group_by(perm_type, model_name,  chr,  in_lab, variable) %>%
  summarize(mean_rnvp = mean(median.val),
            mean_95  = mean(upper.95),
            mean_05  = mean(lower.05)) ->
  sum_dat_for_plot

sum_dat_for_plot %>%
  ggplot(aes(
    x=factor(chr, levels = c("2L", "2R", "3L", "3R")),
    shape = perm_type,
    y=-log10(mean_rnvp),
    ymin=-log10(mean_05),
    ymax=-log10(mean_95),
    color=in_lab,
    fill = in_lab
  )) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.7)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(22, 21)) +
  geom_point(size = 3.0,  color = "black", position = position_dodge(width = 0.7)) +
  theme_bw() +
  xlab("chr") +
  theme(legend.position = "top") +
  facet_wrap(~factor(variable, levels = c("rnp.binom.p", "wZa.p")), ncol = 2) -> perm_summ_plot

ggsave(perm_summ_plot, file = "perm_summ_plot.pdf", w = 6, h = 2.7)




### Plot the zoom in to chr 2 region

all.mod.out %>%
  filter(chr == "2L",
         #perm == 0, 
         label %in% c(
                      #"F","D",
                      "E") ) %>%
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
  #facet_grid(model_name~.) +
  theme_bw() +
  theme(legend.position = "none") ->
  model_E_plot

ggsave(model_E_plot, file = "model_E_plot.pdf", h = 3, w = 6)


############# plot both

ggsave(perm_summ_plot / model_E_plot, file = "summod_model_E_plot.pdf", h = 6, w = 6)


############

#############
load("../4.GML_plots/joint.analysis.perm.real.processed.Rdata")
joint.analysis.perm.real.processed.2L = filter(joint.analysis.perm.real.processed, chr.x == "2L")


##find peaks
joint.analysis.perm.real.processed.2L %>%
  filter(rnp.binom.p < quant_h, invName == "2Lt",
         abs(log10(rnp.binom.p)) > abs(log10(1e-25)) & abs(log10(rnp.wZa.p)) > abs(log10(1e-25)),
         ) %>%
  arrange(start) 

###Peaks are
#4706170-4755827
#5155634-5205808
#6155704-6205728-6255610-6305685-6355769
#6805780
#9555637-9606002

min_peak = c(4706170-1e5, 5155634-0.5e5, 6155704, 6805780, 9555637 )
max_peak = c(4755827, 5205808, 6355769, 6805780, 9606002 )

data.frame(Win.start=min_peak, Win.stop=max_peak) %>%
  mutate(mid.point = (Win.start+Win.stop)/2,
         chr = "2L",
         inv = "2Lt") %>%
  mutate(start=mid.point-1e5,
         stop=mid.point+1e5) %>%
  mutate(lenght = stop-start) %>%
  mutate(WinName = paste(RoundTo(mid.point, 1e5, "floor")/1e6 , "Mb", sep = " " ) ) ->
  PEAKS_for_ANALYSIS

save(PEAKS_for_ANALYSIS, file = "PEAKS_for_ANALYSIS.Rdata")
load("./PEAKS_for_ANALYSIS.Rdata")
#############
ggplot() +
  geom_rect(data = PEAKS_for_ANALYSIS, 
            aes(xmin=start,
                xmax=stop,
                ymin=-65,
                ymax=50,
                fill = WinName), 
            #linetype="solid", color = "purple", 
            alpha = 0.4, fill = "gold1" 
  ) +
  geom_hline(yintercept  = -log10(1e-25), linetype = "dashed", size = 0.45) +
  geom_hline(yintercept  = log10(1e-25), linetype = "dashed", size = 0.45) +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=start), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], 
             aes(xintercept=stop), linetype="solid") +
  geom_line(data = joint.analysis.perm.real.processed.2L, 
            aes(x=start, y= -log10(rnp.binom.p)), 
            color = "red", alpha = 0.9) +
  geom_ribbon(data = joint.analysis.perm.real.processed.2L, 
              aes(x=start, ymin= log10(rnp.wZa.p)), ymax = 0,
              fill = "blue", alpha = 0.4) +
  geom_ribbon(data = joint.analysis.perm.real.processed.2L, 
              aes(x=start, ymax= -log10(quant_h), ymin = -log10(quant_l),), 
              fill = "purple" , alpha = 0.4) +
  #geom_point(data = joint.analysis.perm.real.processed.2L, 
  #           aes(x=start, y= -log10(rnp.binom.p)), 
  #           color = "orange", alpha = 0.3, size = 1.1) +
  geom_point(data = filter(joint.analysis.perm.real.processed.2L, rnp.binom.p < quant_h ),  
             aes(x=start, y= -log10(rnp.binom.p)), 
             alpha = 0.9, shape = 23, size = 0.7) +
  theme_bw() +
  xlim(0,22605762) +
  ylim(-65,50) +
  facet_wrap(~chr.x, scales = "free_x") ->
  plot_perm_real_final_with_wZA

ggsave(plot_perm_real_final_with_wZA, file = "plot_perm_real_final_with_wZA.pdf", h = 4, w = 8)

####### Tajimas D
####### Tajimas D
####### Tajimas D
####### Tajimas D
####### Tajimas D

## old data
## old data
## old data
## old data
load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/0.OLD_ANALYSIS.preAprl1/pi_D_datforplot.Rdata")

D_dat = sub_pi_d_parsed %>%
  filter(resolution == "W_100000",
         variable == "TajimaD") %>% 
  dplyr::select(BIN_START, pop, type, variable, value)

load("./PEAKS_for_ANALYSIS.Rdata")
PEAKS_for_ANALYSIS

D_dat %>%
  filter(type != "all", pop == "cm") -> Ddat_plot

  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_line(
      data=Ddat_plot
      ,aes(
      x=BIN_START,
      y=value,
      color=type,
    )) + 
  geom_rect(data = PEAKS_for_ANALYSIS, 
            aes(xmin=start,
                xmax=stop,
                ymin=-2.2,
                ymax=2,
                fill = WinName), 
            #linetype="solid", color = "purple", 
            alpha = 0.4, fill = "gold1" 
            ) +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="solid") +
  theme_bw() +
  xlim(0,22605762) +
  ylab("Tajima's D") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  old_D_plot

ggsave(old_D_plot, file = "old_D_plot.pdf", w = 7, h = 2.4)



####### LD
####### LD
####### LD
####### LD
####### LD

#### LD plot
load("../7.LD/merged.ld.Rdata")

#ld_df %>%
#  group_by(pair_id) %>%
#  mutate(min_pos = min(BP_A, BP_B)) %>%
#  filter(SNP_A != SNP_B ) %>%
#  group_by(min_pos) %>%
#  summarize(Median_R2 = median(R2),
#            Mean_R2 = mean(R2)) ->
#  summarize_ld_df
#
#save(mean_ld_df, file = "mean_ld_df.forPlot.Rdata")

ld_df %>%
  filter(R2 > 0.6) %>% 
  group_by(pair_id) %>%
  mutate(min_pos = min(BP_A, BP_B)) %>% 
  mutate(anchor_snp = case_when(min_pos == BP_A ~  BP_A,
                                min_pos != BP_A ~  BP_B)) %>%
  filter(SNP_A != SNP_B ) %>%
  group_by(anchor_snp) %>%
  summarize(N= n()) ->
  summarize_N_df


load("../7.LD/mean_ld_df.forPlot.Rdata")
mean_ld_df %>% 
  filter(min_pos < 13e6,
         min_pos > 2.9e6) %>%
  mutate(Pos_round = RoundTo(min_pos, 50000)) %>%
  group_by(Pos_round) %>%
  summarize(Win_R2 = mean(Mean_R2),
            Win_R2sd = sd(Mean_R2)) -> summ_dat

  ggplot() +
    geom_rect(data = PEAKS_for_ANALYSIS, 
              aes(xmin=start,
                  xmax=stop,
                  ymin=0,
                  ymax=0.25,
                  fill = WinName), 
              #linetype="solid", color = "purple", 
              alpha = 0.4, fill = "gold4" 
    ) +
    geom_ribbon(
    data=summ_dat,
    alpha = 0.4,
  aes(
    x=Pos_round,
    y=Win_R2,
    ymin = Win_R2-Win_R2sd,
    ymax = Win_R2+Win_R2sd )) +
    geom_line(data=summ_dat,
              alpha = 0.4
              ,aes(
                x=Pos_round,
                y=Win_R2,
              )) +
    geom_point(data=summ_dat,
              alpha = 0.4
              ,aes(
                x=Pos_round,
                y=Win_R2,
              )) +
    theme_bw() ->  R2_plot_mean

ggsave(R2_plot_mean, file = "R2_plot_mean.pdf", h = 3)

######
ggsave(plot_perm_real_final_with_wZA/old_D_plot/R2_plot_mean, 
       file = "Figure3.new.pdf", w = 6, h = 6)

#### LD among groups -- supp data
#### LD among groups
#### LD among groups
#### LD among groups
#### LD among groups
#### LD among groups
#### LD among groups
setDT(PEAKS_for_ANALYSIS) 
setkey(PEAKS_for_ANALYSIS, start, stop)


ld_df %>%
  mutate(Win_stat_A = case_when(BP_A > 4580998 & BP_A < 4780998 ~ 4.6,
                                BP_A > 5055721 & BP_A < 5255721 ~ 5.1,
                                BP_A > 6155736 & BP_A < 6355736 ~ 6.2,
                                BP_A > 6705780 & BP_A < 6905780 ~ 6.8,
                                BP_A > 9480820 & BP_A < 9680820 ~ 9.5),
         Win_stat_B = case_when(BP_B > 4580998 & BP_B < 4780998 ~ 4.6,
                                BP_B > 5055721 & BP_B < 5255721 ~ 5.1,
                                BP_B > 6155736 & BP_B < 6355736 ~ 6.2,
                                BP_B > 6705780 & BP_B < 6905780 ~ 6.8,
                                BP_B > 9480820 & BP_B < 9680820 ~ 9.5)) ->
  ld_df_winanot

ld_df_winanot$Win_stat_A[is.na(ld_df_winanot$Win_stat_A)] = 0
ld_df_winanot$Win_stat_B[is.na(ld_df_winanot$Win_stat_B)] = 0


ld_df_winanot %>% 
  mutate(winComp = case_when(Win_stat_A == Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A < Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A > Win_stat_B ~ paste(Win_stat_B, Win_stat_A, sep = "_"))) ->
  ld_df_winanot_pair

## Save intermediary file
#save(ld_df_winanot_pair, file = "../7.LD/ld_df_winanot_pair.Rdata")

### Make LD TRIANGLE
### Make LD TRIANGLE
### Make LD TRIANGLE
### Make LD TRIANGLE
### Make LD TRIANGLE
### Make LD TRIANGLE

load("../7.LD/ld_df_winanot_pair.Rdata")
load("./PEAKS_for_ANALYSIS.Rdata")
PEAKS_for_ANALYSIS

### the whole smorgasboard
ggplot() +
  geom_point(data=filter(ld_df_winanot_pair, R2 >= 0.6),
  aes(
    x=BP_A/1e6, y=BP_B/1e6,
    #shape = as.factor(Win_stat_A),
    color=R2, #size = R2
  ),size = 0.09, alpha = 0.1) +
  #+ geom_point( shape = 15, alpha = 0.5) +
  #geom_abline(intercept = 0, slope = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.80) +
  geom_vline(data = PEAKS_for_ANALYSIS, 
             aes(xintercept=mid.point/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.2, size = 3, color = "gold1"
  ) +
  geom_hline(data = PEAKS_for_ANALYSIS, 
             aes(yintercept=mid.point/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.2, size = 3, color = "gold1"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),) ->
  #facet_wrap(~winComp, scale ="free") ->
  ld_triag_smorg

ggsave(ld_triag_smorg, file = "ld_triag_smorg.png", w=3, h=3)



### faceted by windows
ggplot(
  data=filter(ld_df_winanot_pair, R2 >= 0.6, 
              BP_A != BP_B,
              Win_stat_A != Win_stat_B,
              Win_stat_A != 0,Win_stat_B != 0 ),
  aes(
  x=BP_A/1e6, y=BP_B/1e6,
  color=R2, #size = R2
)) +
geom_point(size = 0.3, shape = 15) +
#+ geom_point( shape = 15, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.80) +
  theme_bw() +
  facet_wrap(~winComp, scale ="free") ->
  ld_triag_grid

ggsave(ld_triag_grid, file = "ld_triag_grid.png", w=8, h=7)


#### Make threshold plot
ld_df_winanot_pair %>% .$winComp %>% table %>% names %>% sort -> comb_types
comb_types


PEAKS_for_ANALYSIS %<>%
  separate(WinName, into = c("winComp", "unit"), sep = " ", remove = F)  

#### all snps
ld_df_winanot_pair %>% dim %>% .[1] -> all_pairs
  
ld_df_winanot_pair %>%
  filter(R2 > 0.6) %>%
  summarize(n_r2 = n()) ->
  count_r6_all

#### count r6 by windows
ld_df_winanot_pair %>%
  filter(!winComp %in% comb_types[grep("0",comb_types)]  ) %>% 
  group_by(winComp) %>%
  filter(R2 > 0.6) %>%
  summarize(n_r2 = n()) ->
  count_r26
  
count_r26 %<>%
  separate(winComp, into = c("WinA", "WinB"), sep = "_", remove = F)  

PEAKS_for_ANALYSIS %>%
  dplyr::select(LA = lenght, WinA= winComp) -> PEAKS_A
PEAKS_for_ANALYSIS %>%
  dplyr::select(LB = lenght, WinB= winComp) -> PEAKS_B

count_r6_all/all_pairs -> treshold
count_r26 %>%
  mutate(scaled_r2n = n_r2/2e5) %>%
  ggplot(aes(
    x=winComp,
    y=scaled_r2n
  )) +
  geom_hline(yintercept = treshold$n_r2, size = 0.4, linetype = "dashed") +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() ->
  bp_r2_win

ggsave(bp_r2_win, 
       file = "bp_r2_win.pt.pdf", w = 6, h = 6)


### SUPP DATA
### SUPP DATA
### SUPP DATA
### SUPP DATA
### SUPP DATA

### NEW DATA
pops = c("ME", "PA", "EUROPE", "CM", "AFRICA")
#kars = c("INV", "STD")

root <- "/scratch/yey2sn/Overwintering_ms/16.Haplotypes/"

D_stats = foreach(i=1:length(pops), .combine = "rbind" )%do%{
  message(pops[i])
  
  inv.tmp <- fread( paste(root, "D.", pops[i],  ".W_100000.S_50000.INV.Tajima.D", sep = "" ) )
  inv.tmp %<>% mutate(pop = pops[i], karyo = "INV")
  
  std.tmp <- fread( paste(root, "D.", pops[i],  ".W_100000.S_50000.STD.Tajima.D", sep = "" ) )
  std.tmp %<>% mutate(pop = pops[i], karyo = "STD")
  
  rbind(inv.tmp,std.tmp )
  
}

D_stats %>%
  ggplot(aes(
    x=BIN_START+50000,
    y=TajimaD,
    color=karyo,
    #linetype=pop
  )) +
  geom_hline(yintercept = 0, linetype="solid") +
  #geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  #geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  theme_bw() +
  facet_grid(CHROM~pop, scales = "free_x") +
  geom_line() ->
  d_plot_world

ggsave(d_plot_world, file = "d_plot_world.pdf", h = 6, w = 9)

