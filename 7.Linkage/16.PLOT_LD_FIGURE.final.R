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


###
final.windows.pos = 
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

#### LD plot
load("../7.LD/merged.ld.Rdata")

#### Sneak a peek on the LD of candidates
ld_df %>%
  filter(CHR_A == "2L") %>%
  filter(SNP_A == "2L_5192177_SNP" | SNP_B == "2L_5192177_SNP") ->
  snp.achor.2L_5192177_SNP
inv.maerks = vroom("in2lt_ld_47snps_informative_markers.txt", delim = "_", col_names = F)
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


### count ld pairs
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


ld_df_winanot %>% 
  mutate(winComp = case_when(Win_stat_A == Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A < Win_stat_B ~ paste(Win_stat_A, Win_stat_B, sep = "_"),
                             Win_stat_A > Win_stat_B ~ paste(Win_stat_B, Win_stat_A, sep = "_"))) ->
  ld_df_winanot_pair

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
