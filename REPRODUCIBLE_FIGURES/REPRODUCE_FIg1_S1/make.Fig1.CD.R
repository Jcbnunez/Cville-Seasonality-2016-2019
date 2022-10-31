rm(list = ls())

library(tidyverse)
library(magrittr)
library(forcats)
library(FactoMineR)
library(gtools)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gmodels)
library(MASS)

### panel DC
load("Fig1.panelC.Rdata")

Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  group_by(pop1, bin_date) %>% 
  summarise(FST_mean = mean(FST)) %>%
  dcast(pop1 ~ bin_date) ->
  mean_fst

Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  left_join(mean_fst) %>% 
  ggplot(
    aes(
      x=fct_reorder(pop1, `1.within`),
      y=FST,
      color=bin_date
    )
  ) +
  geom_boxplot(width = 0.7) +
  xlab("Sampling Locale") +
  coord_flip() +
  theme_bw() +
  ylab(expression(italic(F)[ST])) -> 			
  fst_boxplot

ggsave(fst_boxplot,
       file = "fst_boxplot.pdf",
       width = 5,
       height = 4)

### Quantify

Out_comp_vector_samepops %>%
  group_by(bin_date) %>%
  summarise(med = quantile(FST, 0.5 ))

anova(lm(FST ~ bin_date, data = filter(Out_comp_vector_samepops, bin_date %in% c("1.within", "2.Overwinter") )))

### panel D
load("Fig1.panelD.Rdata")

multiyear_samps %>%  
  group_by(pop1, year_diff) %>%
  summarize(Mean = mean(FST),
            SD = sd(FST)) %>%
  ggplot(
    aes(
      x=(year_diff),
      y=(Mean),
      ymin=(Mean)-SD,
      ymax=(Mean)+SD,
      fill=pop1,
      color=pop1,
      #color=as.factor(day_diff_year_scaled)
    ))  + 
  #geom_boxplot(width = 0.4) +
  geom_errorbar(width = 0.1,position=position_dodge(width=0.5)) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_point(shape = 21, size = 2.3, position=position_dodge(width=0.5)) +
  xlab("Number of Winters") +
  theme_bw() +
  xlab(expression(Delta[Years])) +
  ylab(expression(F[ST])) +
  #scale_shape_manual(values = c(21,22, 23)) +
  scale_color_brewer(palette ="Dark2") +
  scale_fill_brewer(palette ="Dark2") ->
  fst_allpop_overwint

ggsave(fst_allpop_overwint,
       file = "fst_allpop_overwint.pdf",
       width = 6,
       height = 2.3)

## Quantify:
anova(lm(FST ~ year_diff*pop1, data = multiyear_samps ))


