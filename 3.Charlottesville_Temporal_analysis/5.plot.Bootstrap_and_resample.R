### Plot resampling
### 
rm(list = ls())

library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
library(data.table)
library(gmodels)
#install_github('tavareshugo/windowscanr')
library(windowscanr)

root <-   "/scratch/yey2sn/Overwintering_ms/3.Cville_PCA/real_resample_out/"
real_in_files = list()
for(i in 1:500){
  tmp <- fread( paste(root, "real.resample.", i, ".txt", sep = "") )
  real_in_files[[i]] = tmp
}
real_df = do.call(rbind, real_in_files)


root <-   "/scratch/yey2sn/Overwintering_ms/3.Cville_PCA/permutation_resample_out/"
perm_in_files = list()
for(i in 1:500){
  tmp <- fread( paste(root, "permuted.resample.", i, ".txt", sep = "") )
  perm_in_files[[i]] = tmp
}
perm_df = do.call(rbind, perm_in_files)

perm_df$chr = "perm"

rbind(real_df,
      perm_df
      ) %>%
  group_by(pc,
           chr,
           snp_n,
           test) %>%
  #summarise(Mean = ci(cor^2)[1],
  #          Low =ci(cor^2)[2],
  #          High = ci(cor^2)[3]) ->
  summarise(Median = median(cor^2),
            IQR25 =quantile(cor^2, 0.25),
            IQR75 = quantile(cor^2, 0.75)) ->
  plot_object

## begin plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
cols = append(cols, "grey")

plot_object %>%
  filter(test == "EFC") %>%
  ggplot(aes(
    x=log10(snp_n),
    y=Median,
    ymin =IQR25,
    ymax =IQR75,
    fill = chr
  )) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  theme_bw() + 
  scale_fill_manual(values = cols) +
  facet_grid(test~pc) ->
  EFC

plot_object %>%
  filter(test == "In(2L)t") %>%
  ggplot(aes(
    x=log10(snp_n),
    y=Median,
    ymin =IQR25,
    ymax =IQR75,
    fill = chr
  )) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  theme_bw() + 
  scale_fill_manual(values = cols) +
  facet_grid(test~pc) ->
  INV

plot_object %>%
  filter(test == "lm_month_2") %>%
  ggplot(aes(
    x=log10(snp_n),
    y=Median,
    ymin =IQR25,
    ymax =IQR75,
    fill = chr
  )) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  theme_bw() + 
  scale_fill_manual(values = cols) +
  facet_grid(test~pc) ->
  lm2m

plot_object %>%
  filter(test == "year") %>%
  ggplot(aes(
    x=log10(snp_n),
    y=Median,
    ymin =IQR25,
    ymax =IQR75,
    fill = chr
  )) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  theme_bw() + 
  scale_fill_manual(values = cols) +
  facet_grid(test~pc) ->
  year

ggsave((year/INV/lm2m/EFC),
       file = "real_data.pdf")
  

  
