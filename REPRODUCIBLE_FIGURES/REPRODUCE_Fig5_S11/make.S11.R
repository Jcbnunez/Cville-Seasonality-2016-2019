### Reproduce fig S11

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(vroom)
library(gmodels)

load("hap.count.data.Rdata")


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

#####
#####
#####
# ------> Investigate Msp300
hap.count.for.plot %>%
  filter(karyo == "std") %>%
  filter(win.st == "win_5.1") %>%
  ggplot(aes(x = mid.bp,
             y= N_UNIQ_HAPS)) +
  geom_line()


