### plot enrichment analyses
### 

### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)

files_o = system("ls *.Rdata", intern = T)

enrich_out = 
foreach(i = 1:4, .combine = "rbind")%do%{
  
  load(files_o[i])
  return(o)
  
}

sets <- data.table(set_id=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


enrich_out %>%
  separate(set, into = c("variable", "set_id"), sep = "_", remove = F ) %>%
  mutate(set_id = as.numeric(set_id)) %>%
  left_join(sets) %>%
  mutate(model_tidy = paste(variable, start, end, sep = "_")) %>%
    ggplot(aes(
      x=chr,
      y=or,
      ymin=lci,
      ymax=uci,
      #color = mod,
      color = inv
    )) +
    geom_hline(yintercept =1) +
    geom_point(size = 2, position=position_dodge(width=0.5)) +
    scale_y_continuous(trans='log2') +
    geom_errorbar(
      width = 0.1,
      position=position_dodge(width=0.5)) +
  theme_bw() +
    facet_grid(model_tidy~mod) -> core_20_plot
  
  ggsave(core_20_plot, file = "core_20_plot.pdf", h = 7)
 
