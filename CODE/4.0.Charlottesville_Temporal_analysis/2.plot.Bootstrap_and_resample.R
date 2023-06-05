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
library(foreach)

###

root <-   "/scratch/yey2sn/old_scra/Overwintering_ms/3.Cville_PCA/Aprl2022_boot_resamp/"
files = system(intern = T, "ls /scratch/yey2sn/old_scra/Overwintering_ms/3.Cville_PCA/Aprl2022_boot_resamp")

perm.input = foreach(i=1:length(files), .combine = "rbind", .errorhandling = "remove")%do%{
  
  message(i/length(files)*100)
  
  fread( 
    paste(root, files[i], sep = "") 
         ) # close fread
  
}

perm.input %>%
  melt(id = c("Pop", "chr", "type",  "jobid", "cov",  "sample_nsnps" )) %>% 
  separate(variable, into = c("variable", "PC"), sep = "_") ->
  perm.input.melt

perm.input.melt %>%
  group_by(Pop,
           chr,
           type,
           cov,
           variable,
           PC) %>%
  summarize(N = n())


perm.input.melt$PC[which(perm.input.melt$type == "Permutation")] = "Perm"

perm.input.melt %>%
  group_by(Pop,
           chr,
           #type,
           sample_nsnps,
           variable,
           PC) %>% 
  summarise(Median = median(value^2),
            IQR25 =quantile(value^2, 0.25),
            IQR75 = quantile(value^2, 0.75),
            IQR05 =quantile(value^2, 0.05),
            IQR95 =quantile(value^2, 0.95),
            Mean = mean(value^2),
            SD = sd(value^2)
            ) ->
  plot_object_dat

save(plot_object_dat, file = "make.FigS4.R")

### plotting all only

plot_object_dat %>% 
  filter( #Pop == "Cha",
         sample_nsnps == 20000) %>%
  ggplot(aes(
    x=variable,
    y=Median,
    ymin=IQR05,
    ymax=IQR95,
    color=PC,
    fill = PC
  )) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.9)) +
  geom_point(size = 1.2, shape = 21, color = "black", position = position_dodge(width=0.9)) +
  theme_bw() +
  coord_flip() +
  facet_grid(Pop~chr) ->
  cvile_plot_resamp

ggsave(cvile_plot_resamp, file = "cvile_plot_resamp.pdf", w=11, h =9)

## plot VA only

plot_object_dat$PC = factor(plot_object_dat$PC)

plot_object_dat %>% 
  filter( Pop == "Cha",
          variable %in% c("Year", "In2Lt") ) %>%
  ggplot(aes(
    x=log10(sample_nsnps),
    y=Median,
    ymin=IQR25,
    ymax=IQR75,
    color=PC,
    fill = PC
  )) +
  geom_ribbon(alpha = 0.1) +
  geom_line(color = "black") +
  #geom_point(size = 1.2, shape = 21, color = "black") +
  theme_bw() +
  facet_grid(chr~factor(variable, levels = c("Year", "In2Lt")  )) ->
  cvile_plot_resamp_va

ggsave(cvile_plot_resamp_va, file = "cvile_plot_resamp_va.pdf", w=4, h =5)
