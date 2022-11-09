#deficiency final figure creation and modeling


library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggrepel)
library(gmodels)
library(lme4)
library(emmeans)
library(patchwork)
registerDoParallel(4)
wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
setwd( wd )

longstartle = readRDS("def.longstartle")
model.data = readRDS("modelingwithgeno.RDS")

#create and save graphs
#start with startle duration and basal activity for each

group.colors <- c(inverted.derived = "blue", standard.ancestral = "red", inverted.ancestral = "steelblue3")

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


#try emmmeans
#need to merge f1.background into geno for this part
#model.data$geno = paste(model.data$geno, model.data$f1.background, sep = "-")

#kno = "k9a"
#make pairwise ~ (categorical factor), var = continues covariate
#model.em = glmer(formula = (zed.obs ~ scaled.minute * geno  + (1| fly.id) + (1|week)), model.data[knockout == kno])

#x = emtrends(model.em, pairwise ~ geno, var = "scaled.minute")#pbkrtest.limit = 6697
#look for emmeans measure of intercepts. 

#emmip(model.em, geno ~ scaled.minute, cov.reduce = range)
#x = x[[1]]
#em.data = as.data.table(x)
#em.data$f1.background = ifelse(grepl("balancer", em.data$geno),"balancer", "deficiency")
#em.data[,newgeno := tstrsplit(geno, "-")[[1]]]
#em.data$knockout = kno
#em.data5 = em.data
#em.total = rbind(em.data1, em.data2, em.data3, em.data4, em.data5)
#saveRDS(em.total, "trends.startle.decay")
em.total = readRDS("trends.startle.decay")

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
ggsave(g3/g1/g2 , file =  "sr.total.partitioned.pdf", width = 9, height = 5)

## Modeling
# use a loop to go through each knockout
knockouts = unique(data.merge$knockout)
out = foreach(f = knockouts) %do% {
 # f = knockouts[1]
  #start with basal act 

mod.ba.full = lmer(formula = (value ~ geno * f1.background + (1 | week)), data = longstartle[phenotypes == "basal.act.min"][knockout == f])

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
