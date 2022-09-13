### Reconstruct boot and resample plot
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

load("make.FigS4.Rdata")


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

#####
write.table(plot_object_dat, 
            file = "TableS4.corelations.all.pops.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")




#### Quantify

plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "Year") %>%
  filter( PC %in% 1) %>% as.data.frame()

plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "Year") %>%
  filter( PC %in% 2) %>% as.data.frame()

plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "Year") %>%
  filter( PC %in% "Perm") %>% as.data.frame()

###
plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "In2Lt") %>%
  filter( PC %in% 1) %>% as.data.frame()

plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "In2Lt") %>%
  filter( PC %in% 2) %>% as.data.frame()

plot_object_dat %>%
  filter(Pop == "Cha") %>%
  filter(sample_nsnps == 20000) %>%
  filter(variable == "In2Lt") %>%
  filter( PC %in% "Perm") %>% as.data.frame()

###
