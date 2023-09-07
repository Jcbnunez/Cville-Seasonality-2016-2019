### plot enrrichment

### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)
library(patchwork)
####

####  create window objects
final.windows.pos = 
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

### load thermal GLM object
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


####### PLot
load("./window.enrich.set.Rdata")
#load("./machado.datasets.enrch.Rdata")

#rbind(mutate(machado.datasets.enrch, set.type = "Machado"  ), 
#      mutate(enrichment.sets, set.type = "DEST")) %>% 
  #filter(analysis_type == "best_model") %>%
  #separate(anchor.model, into = c("model", "resolution.mod", "demo.region"), sep = ";" ) %>%
  enrichment.sets %>%
  mutate(start=win.start,
         end=win.end
  ) %>%      left_join(final.windows.pos) %>%
  ggplot(aes(
    x=win.name,
    y=log2(or),
    ymin=log2(lci),
    ymax=log2(uci),
    color = anchor.model,
    fill = anchor.model,
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(size = 0.5, width = 0.25, position=position_dodge(width=0.5)) +
  geom_point(size = 2.0, shape = 21, position=position_dodge(width=0.5), color = "black") +
  theme_bw() +
  #facet_grid(set.type~., scales = "free_y") +
  theme(legend.pos = "bottom") ->
  enrich.plot

ggsave(enrich.plot, file = "enrich.plot.pdf", w = 5, h = 4)

rbind(mutate(machado.datasets.enrch, set.type = "Machado"  ), 
      mutate(enrichment.sets, set.type = "DEST")) %>% 
  #filter(analysis_type == "best_model") %>%
  #separate(anchor.model, into = c("model", "resolution.mod", "demo.region"), sep = ";" ) %>%
  mutate(start=win.start,
         end=win.end
  ) %>%
  left_join(final.windows.pos) %>%
  ggplot(aes(
    x=win.name,
    y=(st.pr),
    ymin=(st.lci),
    ymax=(st.uci),
    color = anchor.model,
    fill = anchor.model,
  )) +
  geom_errorbar(size = 0.5, width = 0.25, position=position_dodge(width=0.5)) +
  geom_point(size = 2.0, shape = 21, position=position_dodge(width=0.5), color = "black") +
  theme_bw() +
  facet_grid(set.type~., scales = "free_y") +
  theme(legend.pos = "bottom") ->
  dir.plot

ggsave(enrich.plot/dir.plot, file = "dir.plot.pdf", w = 3.5, h = 6)

