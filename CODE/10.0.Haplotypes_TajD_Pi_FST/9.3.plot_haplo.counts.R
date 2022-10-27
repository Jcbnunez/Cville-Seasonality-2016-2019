### Plot haplo counts
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(vroom)
library(gmodels)

inv.haps <- vroom("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/HAPCOUNT.CM.W_100000.S_50000.INV.hapcount")

std.haps <- vroom("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/HAPCOUNT.CM.W_100000.S_50000.STD.hapcount")

inv.haps %>%
  dplyr::select(CHROM, BIN_START,  BIN_END, N_SNP, N_UNIQ_HAPS) %>%
  mutate(karyo = "inv",
  scaled.hap = scale(N_UNIQ_HAPS)) ->
  inv.haps.parsed

std.haps %>%
  dplyr::select(CHROM, BIN_START,  BIN_END, N_SNP, N_UNIQ_HAPS) %>%
  mutate(karyo = "std",
         scaled.hap = scale(N_UNIQ_HAPS)) ->
  std.haps.parsed


#### plot
final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


#After filtering, kept 9 out of 203 Individuals [INV]
#After filtering, kept 152 out of 203 Individuals [STD]

rbind(inv.haps.parsed, std.haps.parsed) %>%
  as.data.frame() %>%
  mutate(N_total = case_when(karyo == "inv" ~ 9*2,
                             karyo == "std" ~ 152*2 ) ) %>%
  mutate(mid.bp = c(BIN_START+BIN_END)/2 ) %>%
  filter(N_UNIQ_HAPS != 0) %>%
  mutate(inv.status = case_when(mid.bp > 2225744 & mid.bp < 13154180 ~ "inv",
                                TRUE ~ "out.inv"),
         win.st = case_when(mid.bp > 2800000 & mid.bp < 3200000 ~ "win_3.1",
                            mid.bp > 4470000 & mid.bp < 4870000 ~ "win_4.7",
                            mid.bp > 4950000 & mid.bp < 5350000 ~ "win_5.1",
                            mid.bp > 6000000 & mid.bp < 6400000 ~ "win_6.1",
                            mid.bp > 6600000 & mid.bp < 7000000 ~ "win_6.8",
                            mid.bp > 9400000 & mid.bp < 9800000 ~ "win_9.6",
                           TRUE ~ "no.win") ) %>%
  mutate(window.label = paste(win.st, inv.status, sep = "_")) -> hap.count.for.plot

### plot
ggplot() +
geom_rect(data = final.windows.pos,
          aes(xmin=start, xmax = end,
              ymin = 0, ymax = 19), 
          alpha = 0.7, fill = "gold") +
  geom_point(data = hap.count.for.plot, aes(
    x=mid.bp,
    y=N_UNIQ_HAPS,
    color = karyo
  )) +
  facet_grid(karyo~. , scales = "free_y") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) +
  theme_bw() +
  scale_color_manual(values = c("firebrick", "steelblue")) ->
  hap.counts

ggsave(hap.counts, file ="hap.counts.pdf", w = 9, h = 3)

#### boxplots
hap.count.for.plot %>%
  ggplot(
    aes(
      x=window.label,
      y=scaled.hap,
      fill = karyo
    )
  ) + 
  geom_boxplot() +
  scale_fill_manual(values = c("firebrick", "steelblue")) ->
  box.plot
ggsave(box.plot, file ="box.plot.pdf",  w = 9, h = 3)

#### Analyses

lm(scaled.hap ~ window.label + karyo, data = hap.count.for.plot) -> w.mod

anova(w.mod)

TukeyHSD(aov(w.mod))$window.label %>%
  as.data.frame() -> TukeyHSD.hapcount

TukeyHSD.hapcount %>% 
  filter(`p adj` < 0.05)

write.table(TukeyHSD.hapcount, 
            file = "TukeyHSD.hapcount.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = T,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")

###
#
save(hap.count.for.plot, file = "hap.count.data.Rdata")


hap.count.for.plot %>%
  mutate(hap.rel = N_UNIQ_HAPS/N_total) %>%
  group_by(karyo) %>%
  summarize(mean = ci(hap.rel)[1],
            lci = ci(hap.rel)[2],
            uci = ci(hap.rel)[3] )


hap.count.for.plot %>%
  mutate(hap.rel = N_UNIQ_HAPS/N_total) %>%
  group_by(karyo, inv.status) %>%
  summarize(mean = ci(hap.rel)[1],
            lci = ci(hap.rel)[2],
            uci = ci(hap.rel)[3] )



