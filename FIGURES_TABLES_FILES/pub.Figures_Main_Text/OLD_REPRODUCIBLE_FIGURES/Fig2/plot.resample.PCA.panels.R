### Plot resample plot


library(tidyverse)
library(magrittr)

load("make.Fig_PCAresample.R")

plot_object_dat$PC = factor(plot_object_dat$PC)

plot_object_dat %>% 
  filter( Pop == "Cha",
          variable %in% c("Year", "In2Lt"),
          sample_nsnps == 1000
          ) %>%
  ggplot(aes(
    x=PC,
    y=Median,
    ymin=IQR25,
    ymax=IQR75,
    fill= variable ,
    shape= variable
    #fill = chr
  )) +
  geom_errorbar(width = 0.1, position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5), size = 2.5) +
  theme_bw() +
  scale_shape_manual(values  = c(21,24)) +
  facet_wrap(chr~., nrow = 1) ->
  cvile_plot_resamp_va_zoom

ggsave(cvile_plot_resamp_va_zoom, 
       file = "zooom_cvile_plot_resamp_va.pdf", w=9, h =2.1)



