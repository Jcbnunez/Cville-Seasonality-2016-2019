####### ====> Plot Genome-wide models
####### 
#### Graph models
#### 
## analysis of models
rm(list = ls())

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(magrittr)
library(ggVennDiagram)
library(patchwork)

results.folder <- "/scratch/yey2sn/old_scra/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis"

files <- system( paste("ls", results.folder, "| grep 'Win'  "), intern = T)

win.all.in =
  foreach(i=1:length(files), 
          .combine = "rbind")%do%{
            print(i)
            load(paste(results.folder, files[i], sep = "/"))
            win.out
          }

win.all.in %>%
  group_by(chr) %>%
  filter(model.pop == "temp.max;2;5.Cville") %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.binom.p, 0.01))  %>%
  mutate(metric = "rnvp") -> rnvp

### Plot
### Plot
### Plot
### Plot
setDT(rnvp)
rnvp %>% group_by(chr) %>% summarize(TN = n()) -> tot.wins
  
rnvp %>% 
  dcast(chr+pos_mean~perm_type, value.var = "uci") %>%
  filter(real < perm) %>%
  group_by(chr) %>% summarize(N = n()) %>%
  left_join(tot.wins) %>%
  mutate(Perc = N/TN)

ggplot() +
  geom_ribbon(
    data=rnvp[perm_type == "perm"],
    aes(
      x=pos_mean/1e6,
      ymax=-log10(uci),
      ymin=0,
    ), alpha = 0.5, color = "grey", fill = "grey", size = 0.1
  ) +
  geom_line(
    data=rnvp[perm_type == "real"],
    aes(
      x=pos_mean/1e6,
      y=-log10(uci),
    ), color = "blue", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  xlab("Genomic Position (Mb)") +
  ylab("Transformed P-value") +
  facet_grid(~chr, scales = "free_x") ->
    rnps.chrs.p3

ggsave(rnps.chrs.p3, file = "rnps.chrs.p3.pdf", w = 7, h = 2.5)



