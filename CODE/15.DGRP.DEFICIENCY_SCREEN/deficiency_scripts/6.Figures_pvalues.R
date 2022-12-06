#Deficiency figure creation and interaction modeling to find significance
#deficiency final figure creation and modeling


library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggrepel)
library(gmodels)
library(patchwork)
registerDoParallel(4)
wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
setwd( wd )

longstartle = readRDS("startle.data.RDS")
em.total = readRDS("emmeans.data.RDS")
data.merge = readRDS('modeling.data.RDS')

#create and save graphs
#start with startle duration and basal activity for each

group.colors <- c(inverted.derived = "red", standard.ancestral = "blue", inverted.ancestral = "green")
g1 = longstartle[phenotypes == "sr.length"] %>%
  
  group_by(geno, f1.background,  phenotypes,knockout,# inversion.st 
  ) %>%
  summarise(mean = ci(value)[1],
            uci = ci(value)[2],
            lci = ci(value)[3]
  ) %>%
  ggplot(aes(
    x=f1.background,
    y=mean,
    ymin=lci,
    ymax=uci,
    color =  geno
  )) +  
  xlab("Gene Deficiency Status") +
  ylab("Startle Duration") +
  scale_color_manual(values = group.colors) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5), show.legend = F) +
  geom_point(position=position_dodge(width = 0.5), show.legend = F) +
  facet_grid(.~knockout, scales = "free_y")+
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

#create graph for basal activity
g2 = 
  longstartle[phenotypes == "basal.act.min"] %>%
  
  group_by(geno, f1.background,  phenotypes,knockout,# inversion.st 
  ) %>%
  summarise(mean = ci(value)[1],
            uci = ci(value)[2],
            lci = ci(value)[3]
  ) %>%
  ggplot(aes(
    x=f1.background,
    y=mean,
    ymin=lci,
    ymax=uci,
    color =  geno
  )) +
  scale_color_manual(values = group.colors) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(.~knockout, scales = "free_y")+
  xlab("Gene Deficiency Status") +
  ylab("Basal Activity") +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
#create graph for modeling

#group.colors <- c(derived = "red", ancestral = "blue")
g3 = 
  ggplot(em.total, aes(x = f1.background, y = scaled.minute.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = newgeno)) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5), show.legend = F) +
  geom_point(position=position_dodge(width = 0.5), show.legend = F) +
  facet_grid(.~knockout, scales = "free_y")+
  theme_bw() + 
  scale_color_manual(values = group.colors) +
  ylab("activity.slope") +
  theme ( 
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

#create patchwork object
ggsave(g3/g1/g2 , file =  "sr.total.partitioned.pdf", width = 14, height = 7)

## Modeling
# use a loop to go through each knockout
knockouts = unique(data.merge$knockout)
out = foreach(f = knockouts) %do% {
   f = knockouts[1]
  #start with basal act 
  
  mod.b4a.full = lmer(formula = (value ~ geno * f1.background + (1 | week)), data = longstartle[phenotypes == "basal.act.min"][knockout == f])
  
  mod.ba.add1 = lmer(formula = (value ~ geno +  f1.background + (1 | week)), data = longstartle[phenotypes == "basal.act.min"][knockout == f])
  
  bastats = anova(mod.ba.full, mod.ba.add1 ,  test="Chisq")
  #find stats for startle length
  mod.srdur.full = lmer(formula = (value ~ geno * f1.background + (1 | week)), data = longstartle[phenotypes == "sr.length"][knockout == f])
  
  mod.srdur.add1 = lmer(formula = (value ~ geno +  f1.background + (1 | week)), data = longstartle[phenotypes == "sr.length"][knockout == f])
  
  sr.dur.stats = anova(mod.srdur.full, mod.srdur.add1,  test="Chisq")
  #find stats for slope decay
  mod.slope.full = lmer(formula = (zed.obs ~ scaled.minute * geno * f1.background  + (1| fly.id) + (1|week)), data.merge[knockout == f])
  
  mod.slope.add1 = lmer(formula = (zed.obs ~ scaled.minute + geno * f1.background  + (1| fly.id) + (1|week)), data.merge[knockout == f])
  
  slope.stats = anova(mod.slope.full, mod.slope.add1,  test="Chisq")
  o = data.frame(
    knockout = f,
    basal.act.df = unlist(bastats[2,7]),
    basal.act.chisq = unlist(bastats[2,6]),
    basal.act.pval = unlist(bastats[2,8]),
    sr.dur.df = unlist(sr.dur.stats[2,7]),
    sr.dur.chisq = unlist(sr.dur.stats[2,6]),
    sr.dur.pval = unlist(sr.dur.stats[2,8]),
    slope.df = unlist(slope.stats[2,7]),
    slope.chisq = unlist(slope.stats[2,6]),
    slope.pval = unlist(slope.stats[2,8])
  )
}
stat.out = rbindlist(out)
write.csv(stat.out, file = "deficiency.modelingstats.csv")
saveRDS(stat.out, "stat.summary.sr")
