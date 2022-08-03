### Panels 3B and 3C
### 
rm(list = ls())


library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)

### Panel 3B
### 

#save(sub_pi_d_parsed.plot, inv.dt, final.windows.pos, file = "dat.for.3b.Rdata")
load("./dat.for.3b.Rdata")

ggplot() + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="solid") +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -0, ymax = 0.01), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=sub_pi_d_parsed.plot,
    aes(
      x=BIN_START/1e6,
      y=value,
      color = type),
    alpha = 0.9) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,20.5) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") ->
  pi_d_plot_all

ggsave(pi_d_plot_all, file = "pi_d_plot_all.pdf", w = 7, h = 3.5)

### Panel 3C
### 

load("GEVA.2L.win.age.Rdata")

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start, xmax = end,
                ymin = 0, ymax = 190000), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=filter(win.out.geva, V11 != "Dmel_ALL"),
    aes(
      x=mean.pos,
      y=(med.age*0.06666667),
      color = V11
    )) +
  #geom_ribbon() +
  geom_line() +
  ylab("TMRCA gens (x1e6)") +
  xlab("Genomic Position (Mb)") +
  ggtitle("TMRCA years") +
  theme_bw() +
  xlim(0,20.5*1e6) +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  geva.win

ggsave(geva.win, file ="geva.win.pdf", w = 7, h = 3)

