### Collect simultion obs data
### 

library(tidyverse)
library(magrittr)
library(foreach)
library(reshape2)

dat_obs = system("ls obsDat.*", intern = T)

obs_summaries = foreach(i=1:48, .combine = "rbind" )%do%{
  
  load(dat_obs[i])
  data.out
  
}

#2R	15391154	20276334	2RNS


obs_summaries %>%
  as.data.frame() %>%
  filter(win.start > 1e6 & win.start < 12e6 ) %>%
  group_by(Obs_var) %>%
  summarize(#Mean_Val = mean(value),
          Median_Val = median(value),
          ) %>% as.data.frame() ->
  obs_values.final
  
save(obs_values.final, file = "obs_values.final.Rdata")
  
  
  mutate(mid_point = (win.start+win.end)/2) %>%
  filter(Obs_var %in% c("yearDiff=1.within_FST", "yearDiff=2.Overwinter_FST", "yearDiff=3.Multi-Year_FST")  ) %>% 
  ggplot(aes(
    x=mid_point,
    y=as.numeric(value),
    color = Obs_var
  )) +
  geom_line() ->
  test_plot

ggsave(test_plot, file = "test_plot.pdf")


obs_summaries %>%
  as.data.frame() %>%
  mutate(mid_point = (win.start+win.end)/2) %>%
  filter(Obs_var %in% c("DimDAPC=LD1_corr.time.est", "DimDAPC=LD2_corr.time.est")  ) %>% 
  ggplot(aes(
    x=mid_point,
    y=as.numeric(value),
    color = Obs_var
  )) +
  geom_line() ->
  test_plot2

ggsave(test_plot2, file = "test_plot2.pdf")

obs_summaries %>%
  as.data.frame() %>%
  ggplot(aes(
    x=Obs_var,
    y=log10(value)
  )) +
  geom_boxplot() +
  coord_flip() ->
  boxplots.vars

ggsave(boxplots.vars, file = "boxplots.vars.pdf")



