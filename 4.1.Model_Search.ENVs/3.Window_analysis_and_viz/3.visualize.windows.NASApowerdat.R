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

results.folder <- "/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/GLM_omnibus_window_analysis"

files <- system( paste("ls", results.folder, "| grep 'Win'  "), intern = T)

win.all.in =
foreach(i=1:length(files), 
        .combine = "rbind")%do%{
  print(i)
  load(paste(results.folder, files[i], sep = "/"))
  win.out
}

### Plot distributions

win.all.in %>%
  filter(chr == "2L") %>%
  filter(model.pop == "temp.max;2;5.Cville") %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.binom.p, 0.01))  %>%
  mutate(metric = "rnvp") -> rnvp

win.all.in %>%
  filter(chr == "2L") %>%
  filter(model.pop == "temp.max;2;5.Cville") %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.wZa.p, 0.01))  %>%
  mutate(metric = "wza") -> wza

rbind(rnvp, wza ) %>%
  ggplot(
    aes(
      x=pos_mean/1e6,
      y=-log10(uci),
      color=perm_type
    )
  ) +
  geom_line(alpha = 0.7) +
  ylab("P-values") +
  xlab("Genome Position (Mb)") +
  theme_bw() +
  scale_color_manual(values = c("grey","red") ) +
  facet_grid(chr~metric, scales = "free") ->
  cville.temp.rnp.plot

ggsave(cville.temp.rnp.plot, file = "cville.temp.rnp.plot.pdf")

##----> visualize
##

rnvp %>%  filter(perm_type == "real") %>% arrange(uci)

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

#data.frame(
#  start= c(4650065, 5100324, 6100321, 9500286),
#  end= c(4799922, 5349218, 6349489, 9700005)
#) %>% mutate(mid.point = (start/2) + (end/2) ) %>%
#  mutate(win.name = paste("win", round(mid.point/1e6, 1), sep = "."  )) ->
#  final.windows.pos

rbind(rnvp, wza ) %>%
  mutate(uci.mod = case_when(metric == "rnvp" ~ -log10(uci),
                            metric != "rnvp" ~ log10(uci),
                            ) ) -> dat.for.plot

#####
  ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -110, ymax = 110), 
            alpha = 0.7, fill = "gold") +
   geom_line(
     data=dat.for.plot,
     aes(
       x=pos_mean/1e6,
       y=uci.mod,
       color = paste(perm_type, metric)
     ), alpha = 0.5
   ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-115,110) +
  scale_color_manual(values = c("grey50", "grey50", "blue", "purple")) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Genomic Position (Mb)") +
  ylab("Transformed P-value") +
  xlim(0,20.5) +
   geom_vline(xintercept = 2225744/1e6) +
   geom_vline(xintercept = 13154180/1e6) ->
   panel.for.fig3
 
 ggsave(panel.for.fig3, file = "panel.for.fig3.pdf", w = 7, h = 3)

 
 save(final.windows.pos, dat.for.plot, file = "dat.for.panel3A.Rdata")
 ######
 ######
 ######
 

 