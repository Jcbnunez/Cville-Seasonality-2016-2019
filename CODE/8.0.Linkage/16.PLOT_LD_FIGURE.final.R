### PLOT LD FIGURE
### 
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
library(vroom)

win.bp <- 100000
step.bp <- 50000

#### model files
base <- "/project/berglandlab/alan/environmental_ombibus_global"
models = c("temp.max;2;5.Cville")


file <- paste(base, models, paste( models,"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)
out.glm <- get(load(file))
out.glm.2L = out.glm %>% filter(chr == "2L" & perm == 0)

inv.maerks = vroom("in2lt_ld_47snps_informative_markers.txt", 
                   delim = "_", col_names = F)

### load windows
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.2, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

#### LD plot
load("../7.LD/merged.ld.Rdata")
##load("./merged.ld.Rdata")

#### Sneak a peek on the LD of candidates
ld_df %>%
  filter(CHR_A == "2L") %>%
  filter(SNP_A == "2L_5192177_SNP" | SNP_B == "2L_5192177_SNP") ->
  snp.achor.2L_5192177_SNP


inv.maerks %<>% mutate(snp.id = paste(paste(X1, "L", sep = ""),X2,X3, sep = "_"))

snp.achor.2L_5192177_SNP %>%
  filter(SNP_A %in% inv.maerks$snp.id) %>%
  summarize(meanR2 = mean(R2),
            sd = sd(R2))

####
ld_df %>%
  filter(R2 > 0.6 &
          SNP_A !=  SNP_B) -> ld_df_0.6


ggplot() +
  geom_point(data=ld_df_0.6,
             aes(
               x=BP_A/1e6, y=BP_B/1e6,
               #shape = as.factor(Win_stat_A),
               color=R2, #size = R2
             ),size = 0.09, alpha = 0.1) +
  #+ geom_point( shape = 15, alpha = 0.5) +
  #geom_abline(intercept = 0, slope = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "purple", midpoint =  0.80) +
  geom_vline(data = final.windows.pos, 
             aes(xintercept=((start+end)/2)/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.2, size = 3, color = "gold1"
  ) +
  geom_hline(data = final.windows.pos, 
             aes(yintercept=((start+end)/2)/1e6), 
             #linetype="solid", color = "purple", 
             alpha = 0.2, size = 3, color = "gold1"
  ) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank()
        ) ->
  #facet_wrap(~winComp, scale ="free") ->
  ld_triag_smorg

ggsave(ld_triag_smorg, file = "ld_triag_smorg.pdf", w=4, h=3)


#### Characterizing LD landscapes in outlier windows
#### Characterizing LD landscapes in outlier windows
#### Characterizing LD landscapes in outlier windows
#### Characterizing LD landscapes in outlier windows <<---
#### Characterizing LD landscapes in outlier windows
#### Characterizing LD landscapes in outlier windows
#### Characterizing LD landscapes in outlier windows

ld_df %>%
  mutate(Win_stat_A = case_when(BP_A > 2000000 & BP_A < 2400000 ~ "left",
                                BP_A > 12900000 & BP_A < 13300000 ~ "right",
                                BP_A > 2800000 & BP_A < 3200000 ~ "3.1",
                                BP_A > 4470000 & BP_A < 4870000 ~ "4.7",
                                BP_A > 4920000 & BP_A < 5320000 ~ "5.1",
                                BP_A > 6000000 & BP_A < 6400000 ~ "6.1",
                                BP_A > 6600000 & BP_A < 7000000 ~ "6.8",
                                BP_A > 9400000 & BP_A < 9800000 ~ "9.6"),
         Win_stat_B = case_when(
                                BP_B > 2000000 & BP_B < 2400000 ~ "left",
                                BP_B > 12900000 & BP_B < 13300000 ~ "right",
                                BP_B > 2800000 & BP_B < 3200000 ~ "3.1",
                                BP_B > 4470000 & BP_B < 4870000 ~ "4.7",
                                BP_B > 4920000 & BP_B < 5320000 ~ "5.1",
                                BP_B > 6000000 & BP_B < 6400000 ~ "6.1",
                                BP_B > 6600000 & BP_B < 7000000 ~ "6.8",
                                BP_B > 9400000 & BP_B < 9800000 ~ "9.6")
                                ) ->
  ld_df_winanot

ld_df_winanot$Win_stat_A[is.na(ld_df_winanot$Win_stat_A)] = 0
ld_df_winanot$Win_stat_B[is.na(ld_df_winanot$Win_stat_B)] = 0

###
ld_df_winanot %>%
  filter(R2 > 0.6) %>%
  group_by(pos=BP_A, Win_stat_A) %>%
  summarize(N = n()) %>%
  group_by(Win_stat_A) %>%
  slice_max(N, n=3) %>%
  left_join( dplyr::select(out.glm.2L, pos, outP= p_lrt) ) %>%
  filter(!is.na(outP)) %>%
  slice_min(outP) -> top_hits

ld_df_winanot %>%
filter(BP_A %in% top_hits$pos & BP_A %in% top_hits$pos) %>%
  left_join( dplyr::select(out.glm.2L, BP_B=pos, P_B= p_lrt)) %>%
  mutate(p_lrt.tresh = case_when(
    P_B < 0.1 & P_B >= 0.01  ~ "p-1",
    P_B < 0.01 & P_B >= 0.001 ~ "p-2",
    P_B < 0.001 ~ "p-3",
    P_B > 0.1 ~ "NS")) %>%
  filter(!is.na(p_lrt.tresh)) %>%
  group_by(Win_stat_A, Win_stat_B,p_lrt.tresh) %>%
  summarize(mR2 = mean(R2)) %>%
  filter(Win_stat_B %in% c("3.1","4.7","5.1","6.1","6.8","9.6")) %>% 
  filter(Win_stat_A %in% c("3.1","4.7","5.1","6.1","6.8","9.6")) %>%
  filter(Win_stat_A != Win_stat_B ) %>%
  group_by(p_lrt.tresh) %>%
  summarize(mR2 = mean(mR2))
  


###
#anch = "6.1"
#foreach(anch = c("left","right","3.1","4.7","5.1","6.1","6.8","9.6"))%do%{
ld_df_winanot %>%
  filter(R2 > 0.2) %>% 
  mutate(anchor = anch ) %>%
  mutate(pos_low = ifelse(BP_A > BP_B, BP_B, BP_A ) ) %>%
  mutate(pos_hig = ifelse(BP_A > BP_B, BP_A, BP_B ) ) ->
  window.ex

wondow.o = window.ex %>% filter(R2 > 0.7) %>%
  filter(Win_stat_A != 0 | Win_stat_B != 0) %>%
  filter(BP_A %in% top_hits$pos ) 

  ggplot() +
    geom_jitter(
    data = window.ex,
    aes(
      x=pos_low,
      y=pos_hig,
      color = R2
    ), size = 0.5, shape = 20
  ) + geom_jitter(
    data = wondow.o,
    aes(
      x=pos_low,
      y=pos_hig,
      fill = R2,
      color = R2
    ), shape = 21, size = 2.2
  ) +
  scale_color_gradient2(low = "blue", 
                        high = "red",
                        midpoint = 0.5) +
  scale_fill_gradient2(low = "blue", 
                          high = "red",
                       midpoint = 0.5) +
    theme_classic() ->
  test_triag
ggsave(test_triag, file ="test_triag.png")


outlier.ld =
foreach(anch = top_hits$pos, .combine = "rbind")%do%{
    
ld_df_winanot %>%
  filter(BP_A == anch ) %>%
  mutate(anchor = anch ) %>%
  mutate(pos =  BP_B ) %>%
  left_join(out.glm.2L[,c("pos","p_lrt")]) %>%
  mutate(outlier = ifelse( BP_A == pos, BP_B, BP_A  ) ) %>%
  left_join( dplyr::select(out.glm.2L, outlier=pos, outP= p_lrt) ) %>%
  filter(outP == min(outP, na.rm = T)) %>%
  mutate(p_lrt.tresh = case_when(
  p_lrt < 0.1 & p_lrt >= 0.01  ~ "p-1",
  p_lrt < 0.01 & p_lrt >= 0.001 ~ "p-2",
  p_lrt < 0.001 ~ "p-3"
  #& p_lrt >= 0.0001 ~ "p-3",
  #p_lrt < 0.0001 ~ "p-4",
  #p_lrt < 0.0001 & p_lrt >= 1e-5 ~ "p-4",
  #p_lrt < 1e-5 ~ "p-5"
  #  p_lrt < 0.05 ~ "S",
    TRUE ~ "NS"
  )) -> r2_tmp
  
  wins =
    data.table(#win=win.i,
      start=seq(from=min(r2_tmp$BP_B), 
                to=max(r2_tmp$BP_B)-win.bp, by=step.bp),
      end=seq(from=min(r2_tmp$BP_B), 
              to=max(r2_tmp$BP_B)-win.bp, by=step.bp) + win.bp)
  
  setDT(wins)
  
  win.ld.out <- foreach(win.i=1:dim(wins)[1], 
                        .errorhandling = "remove",
                        .combine = "rbind"
  )%do%{
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    win.tmp <- r2_tmp %>%
      filter(BP_B > wins$start[win.i] & BP_B < wins$end[win.i] ) 
    
    win.tmp %>%
      group_by(p_lrt.tresh, anchor, Win_stat_A) %>%
      summarize(R2_m = mean(R2),
                R.90 = sum(R2 >= 0.9),
                R.60 = sum(R2 >= 0.6),
      ) %>%
      mutate( start=wins$start[win.i],
              end=wins$end[win.i]
      )
    
  }
  
}
save(outlier.ld, file = "outlier.ld.Rdata")
load("outlier.ld.Rdata")

outlier.ld %>%
filter(p_lrt.tresh %in% c("p-3","NS")) %>%
filter(!Win_stat_A %in% c("left","right",0) ) -> outlier.ld.plot

ggplot()  +
  theme_bw() +
  geom_vline(data = final.windows.pos, 
             aes(xintercept=((start+end)/2)), 
             #linetype="solid", color = "purple", 
             alpha = 0.5, size = 3, color = "gold1") + 
geom_smooth(data = outlier.ld.plot,
            aes(
              x=(start+end)/2,
              y=R2_m,
              linetype = p_lrt.tresh,
              color = Win_stat_A
            ),alpha = 0.1) -> anchor_plot
ggsave(anchor_plot, file = "anchor_plot.pdf", w = 8, h = 2.3)

###
###
###
###
###


ld_df_winanot %>% 
  mutate(winComp = case_when(Win_stat_A == Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A < Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A > Win_stat_B ~ paste(Win_stat_B, Win_stat_A, sep = "_"))) ->
  ld_df_winanot_pair

ld_df_winanot_pair$winComp %>% unique() %>% grep("[right,left]", .) -> bkpts
ld_df_winanot_pair$winComp %>% unique() %>% .[bkpts] -> labels

setDT(ld_df_winanot_pair)

inv2ltsvms = inv.maerks$X2

ld_df_winanot_pair %>%
  mutate(marker = case_when(BP_A %in% inv2ltsvms | BP_B %in% inv2ltsvms ~ "svm",
                            TRUE ~ "noSvm") ) ->
  inv2ltsvms.svm
setDT(inv2ltsvms.svm)

inv2ltsvms.svm[winComp == "left_right"][R2 > 0.97] -> inversion.markers

unique(inversion.markers$BP_A, inversion.markers$BP_B) -> high.ld.inv.markers

ld_df_winanot_pair %>%
  filter(BP_A %in% high.ld.inv.markers | BP_B %in% high.ld.inv.markers ) ->
  markers_to_inv

ld_df_winanot_pair %>%
  filter(Win_stat_A %in% c("0","3.1","4.7","5.1","6.1","6.8","9.6") &
          Win_stat_B %in% c("0","3.1","4.7","5.1","6.1","6.8","9.6")) ->
    markers_to_window

####
rbind(markers_to_inv, markers_to_window) -> LD_filtered_dt

### make the landscape of LD across the inversion
LD_filtered_dt %>%
  mutate(delta_bp = abs(BP_A-BP_B)) ->
  LD_filtered_dt

LD_filtered_dt %>%
  group_by(Win_stat_A,delta_bp) %>%
  summarize(meanr2 = mean(R2, na.rm = T)) ->
  summarized_r2

setDT(summarized_r2)

save(summarized_r2, file = "summarized_r2.Rdata")
load("summarized_r2.Rdata")



## prepare windows
## 
wins =
data.table(#win=win.i,
           start=seq(from=min(summarized_r2$delta_bp), 
                     to=max(summarized_r2$delta_bp)-win.bp, by=step.bp),
           end=seq(from=min(summarized_r2$delta_bp), 
                   to=max(summarized_r2$delta_bp)-win.bp, by=step.bp) + win.bp)

wins = rbind(data.frame(#win=1, 
                        start = 0, end = 10), 
             wins, 
             data.frame(#win=1, 
                        start = 11142064-10, end = 11142064)
             )
setDT(wins)


win.ld.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  win.tmp <- summarized_r2 %>%
    filter(delta_bp > wins$start[win.i] & delta_bp < wins$end[win.i] ) 
  
  win.tmp %>%
    group_by(Win_stat_A) %>%
    summarize(R2_m = mean(meanr2),
              R.90 = sum(meanr2 >= 0.9),
              R.60 = sum(meanr2 >= 0.6),
              ) %>%
    mutate( start=wins$start[win.i],
            end=wins$end[win.i]
            )
}
save(win.ld.out, file = "win.ld.summary.out.Rdata")

load("win.ld.summary.out.Rdata")

win.ld.out %>%
  mutate(type = case_when(Win_stat_A %in% c(0,3.1,4.7,5.1,6.1,6.8,9.6) ~ "Window",
                          Win_stat_A %in% c("left","right") ~ "Breakpoint")) %>%
  ggplot(aes(
    x=(start+end)/2,
    y=R2_m,
    color = Win_stat_A
  )) +
  geom_line() +
  facet_grid(type~.)->
  r2plotwin

ggsave(r2plotwin, file = "r2plotwin.pdf", w = 6, h = 8)


win.ld.out %>%
  mutate(type = case_when(Win_stat_A %in% c(0,3.1,4.7,5.1,6.1,6.8,9.6) ~ "Window",
                          Win_stat_A %in% c("left","right") ~ "Breakpoint")) %>%
  filter(type == "Breakpoint") %>%
  ggplot(aes(
    x=(start+end)/2,
    y=R2_m,
    color = Win_stat_A
  )) +
  geom_line() +
  facet_grid(Win_stat_A~type)->
  breakpoints

ggsave(breakpoints, file = "breakpoints.pdf", w = 6, h = 8)




### compare the outlier SNPs to breakpoints -- use the filtered DT
### compare the outlier SNPs to breakpoints
### compare the outlier SNPs to breakpoints
### compare the outlier SNPs to breakpoints

ld_df_winanot_pair %>%
  filter(winComp %in% labels) %>%
  filter(!winComp %in% c("0_left","0_right","right_right","left_left","left_right") ) %>%
  mutate(pos = case_when(Win_stat_A %in% c("right","left") ~  BP_B,
                         Win_stat_B %in% c("right","left") ~  BP_A,))->
  ld_df_winanot_pair.flt

left_join(ld_df_winanot_pair.flt, out.glm.2L) -> ld_df_winanot_pair.flt.glm
setDT(ld_df_winanot_pair.flt.glm)

ld_df_winanot_pair.flt.glm %<>%
  separate(winComp, remove = F, into = c("outloc","bkpt"), sep = "_" )

ld_df_winanot_pair.flt.glm %>%
  mutate(p_lrt.tresh = case_when(
    p_lrt < 0.1 & p_lrt >= 0.01  ~ "p-1",
    p_lrt < 0.01 & p_lrt >= 0.001 ~ "p-2",
    p_lrt < 0.001 & p_lrt >= 0.0001 ~ "p-3",
    p_lrt < 0.0001 & p_lrt >= 1e-5 ~ "p-4",
    p_lrt < 1e-5 ~ "p-5",
  )) ->
  ld_df_winanot_pair.flt.glm.scored

save(ld_df_winanot_pair.flt.glm.scored, file = "ld_df_winanot_pair.flt.glm.scored.Rdata")

ld_df_winanot_pair.flt.glm.scored %>%
  filter(!is.na(p_lrt.tresh)) %>%
  ggplot(aes(
    x=outloc,
    y=R2,
    fill=p_lrt.tresh
  )) +
  geom_boxplot(outlier.shape = NA)  ->
  rnp.tresh.plot

ggsave(rnp.tresh.plot, file = "rnp.tresh.plot.pdf", w = 7, h = 4 )

##### EXTRA ANALYSES
##### EXTRA ANALYSES
##### EXTRA ANALYSES
##### EXTRA ANALYSES
##### EXTRA ANALYSES
##### EXTRA ANALYSES
##### EXTRA ANALYSES

##### Counting LD pairs
ld_df_winanot_pair %>% .$winComp %>% table %>% names %>% sort -> comb_types
comb_types

ld_df_winanot_pair %>%
  filter(!winComp %in% comb_types[grep("0",comb_types)]  ) %>% 
  group_by(winComp) %>%
  filter(R2 > 0.6) %>%
  summarize(n_r2 = n()) ->
  count_r26

count_r26 %<>%
  separate(winComp, into = c("WinA", "WinB"), sep = "_", remove = F)  

#PEAKS_for_ANALYSIS %>%
#  dplyr::select(LA = lenght, WinA= winComp) -> PEAKS_A
#PEAKS_for_ANALYSIS %>%
#  dplyr::select(LB = lenght, WinB= winComp) -> PEAKS_B

#### all snps
ld_df_winanot_pair %>% dim %>% .[1] -> all_pairs

ld_df_winanot_pair %>%
  filter(R2 > 0.6) %>%
  summarize(n_r2 = n()) ->
  count_r6_all

count_r6_all/all_pairs -> treshold

count_r26 %>%
  filter(WinA != WinB) %>%
  filter(!winComp %in% c("left_right", "right_left")) %>%
  mutate(scaled_r2n = n_r2/2e5) %>%
  ggplot(aes(
    x=fct_reorder(winComp, scaled_r2n),
    y=scaled_r2n
  )) +
  geom_hline(yintercept = treshold$n_r2, size = 0.4, linetype = "dashed") +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() ->
  bp_r2_win

ggsave(bp_r2_win, 
       file = "bp_r2_win.pt.pdf", w = 6, h = 6)
