library(tidyverse)
library(magrittr)
library(data.table)
library(gmodels)
library(vroom)


inv.tb <- vroom("InversionsMap_hglft_v6_inv_startStop.txt")
names(inv.tb)[1] = "chr"

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6,
         chr = "2L")

ggplot() +
  #geom_rect(data=outlier_haplowins, aes(xmin=start/1e6, xmax=end/1e6, ymin=-1, ymax=100),fill ="lightgoldenrod1", alpha=.5) +
  geom_rect(data = inv.tb,
            aes(xmin=start/1e6, xmax = stop/1e6,
                ymin = 0, ymax = 35, fill = invName), 
            alpha = 0.7) +
  #geom_line(data = dt.window[perm.st == "permutation"],
  #           aes(x=I(start/2+end/2)/1e6, y=avg.N), color="black") +
  #geom_line(data = dt.window[perm.st == "permutation"],
  #           aes(x=I(start/2+end/2)/1e6, y=lowerbound), color="black", linetype = "dashed") +
  #geom_vline(xintercept = 2225744/1e6) +
  #geom_vline(xintercept = 13154180/1e6) +
  geom_ribbon(data = dt.window[perm.st == "permutation"],
              aes(x=I(start/2+end/2)/1e6, ymax= upperbound, ymin = 0), fill="grey", linetype = "solid", alpha = 0.5) +
  geom_point(data=dt.window[perm.st == "observed"],
             aes(x=I(start/2+end/2)/1e6, y=avg.N, color = sig.v.per), #color="firebrick4",
             size = 2) +
  theme_bw()+
  xlab("Genome Position (Mb)") +
  ylab("Phenotypes (N)") +
  facet_grid(~chr, scales = "free_x")
  theme(legend.position = "none") +
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos") 
  
  
  -> PanelC

ggsave(PanelC, file =  "./new.model2.5.fdr.pdf", w = 7, h = 3)
