### Code to process ld data
### 

library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)
#parallel computing in R
library(foreach)
library(doMC)
library(car)
library(DescTools)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem
###
##########

#### load
#save(ld_list, file = "INV.merged.ld.Rdata")
load("./INV.merged.ld.Rdata")

ld_list %>%
  filter(R2 > 0.5,
         BP_A < BP_B)  %>% 
ggplot(
  aes(
    x=BP_A,
    y=BP_B,
    color=R2
  )
) + geom_point(size = 0.09, alpha = 0.1, shape = 15) +
  theme_classic() +
  scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.80) ->
  inv_ld_plot

ggsave(inv_ld_plot, file = "inv_ld_plot.png")

##########
##########
##########
##########
##########

ld_list %>%
  mutate(#dist_bin = RoundTo(BP_diff, 100000, "floor"),
         BPA_binned = RoundTo(BP_A, 100000, "floor"),
  ) -> 
  ld_dat_merged_binned

###
ld_tresh_count = foreach(i=seq(from=0.5, to = 0.9, by = 0.1), .combine = "rbind")%do%{
  
  message(i)
  
ld_dat_merged_binned %>%
  group_by(BPA_binned) %>%
  filter(R2 > i) %>%
  summarize(N_i = n()) ->
  tmp_obj_flt

ld_dat_merged_binned %>%
  group_by(BPA_binned) %>%
  summarize(N_all = n()) ->
  tmp_obj_all

left_join(tmp_obj_flt, tmp_obj_all) %>%
  mutate(R2_tresh = i,
         N_freq = N_i/N_all ) ->
  ld_dat_merged_binned_mean_per_pos

  return(ld_dat_merged_binned_mean_per_pos)
  
}


####
ld_tresh_count %>%
  ggplot(
    aes(x=BPA_binned,
        y=N_i,
        color=as.factor(R2_tresh)
        )
  ) + 
  geom_line( stat = "identity") +
  xlab("Window Start") +
  facet_wrap(~R2_tresh, ncol = 1, scales = "free") ->
  ld_plots_means

ggsave(ld_plots_means, file = "INV.ld_plots_means.pdf")





