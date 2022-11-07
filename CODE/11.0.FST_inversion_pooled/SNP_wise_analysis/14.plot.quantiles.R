### Plot quantile set
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doParallel)
library(viridis)
#####
in2Lt.markers <- vroom("./inv2L_correlated_markers_Dm3.txt")
quantile.collect.folder = "/scratch/yey2sn/old_scra/Overwintering_ms/19.inv.fst/quantile.ind.snps.out"

###

files.vec.inv = system(paste("ls ",  quantile.collect.folder),
                   intern= T)

collect.quant.inv =
  foreach(i = 1:length(files.vec.inv),
          .combine = "rbind")%do%{
            
            message(paste(i, length(files.vec.inv), sep = " | "))
            tmp <- get(load(paste(quantile.collect.folder, files.vec.inv[i], sep = "/" )))
            
          }

setDT(collect.quant.inv)

collect.quant.inv %<>%
  mutate(SNP_id = anchor.snp)

left_join(collect.quant.inv, in2Lt.markers, by = "SNP_id") %>%
  mutate(type = case_when(cor_rank >= 0.75 & cor_rank < 0.77 ~ "r0.75-0.77",
                          cor_rank >= 0.77 & cor_rank < 0.80 ~ "r0.77-0.80",
                          cor_rank >= 0.80 & cor_rank < 0.83 ~ "r0.80-0.83",
                          cor_rank >= 0.83 & cor_rank < 0.85 ~ "r0.83-0.85",
                          cor_rank >= 0.85 & cor_rank < 0.87 ~ "r0.85-0.87",
                          cor_rank >= 0.87 & cor_rank < 0.90 ~ "r0.87-0.90",
                          cor_rank >= 0.90 & cor_rank < 0.93 ~ "r0.90-0.93",
                          cor_rank >= 0.93 & cor_rank < 0.95 ~ "r0.93-0.95",
                          cor_rank >= 0.95 & cor_rank < 0.97 ~ "r0.95-0.97",
                          cor_rank >= 0.97 & cor_rank < 1.00 ~ "r0.97-1.00",
                          TRUE ~ "r<0.75")) ->
  collect.quant.inv.annot

collect.quant.inv.annot %<>%
  mutate(analysis.set = case_when(year_diff %in% 0:2 ~ "time",
                                  TRUE ~ "space"))


#####
#####
##### ---> GLM set
#####
#####
#####
glm.snps <- get(load("./Cville.GLM.rnp5.snps.Rdata"))
quantile.glm.collect.folder = "/scratch/yey2sn/old_scra/Overwintering_ms/19.inv.fst/quantile.ind.snps.out.glm"

files.vec.glm = system(paste("ls ",  quantile.glm.collect.folder),
                       intern= T)

collect.quant.glm =
  foreach(i = 1:length(files.vec.glm),
          .combine = "rbind")%do%{
            
            message(paste(i, length(files.vec.glm), sep = " | "))
            tmp <- get(load(paste(quantile.glm.collect.folder, files.vec.glm[i], sep = "/" )))
            ##tmp.2 <- get(load("2L_2800243_SNP.quantile.Rdata"))
          }
setDT(collect.quant.glm)


collect.quant.glm %<>%
  mutate(SNP_id = anchor.snp) %>%
  separate(anchor.snp, remove = F, into = c("chr", "pos", "var.type"), sep = "_")

  collect.quant.glm$pos = as.numeric(collect.quant.glm$pos)  
  
left_join(collect.quant.glm, glm.snps) %>% 
  mutate( type =case_when(
    pos > 2800000 & pos < 3200000 ~ "glm.w3.1",
    pos > 4470000 & pos < 4870000 ~ "glm.w4.7",
    pos > 4920000 & pos < 5320000 ~ "glm.w5.1",
    pos > 6000000 & pos < 6400000 ~ "glm.w6.1",
    pos > 6600000 & pos < 7000000 ~ "glm.w6.8",
    pos > 9400000 & pos < 9800000 ~ "glm.w9.6",
    TRUE ~ "noWin")) %>%
  mutate(rnp.thresh = 
    case_when(
    rnp <= 0.05 & rnp > 0.01 ~ "5%",
    rnp <= 0.01 & rnp > 0.001 ~ "1%",
    rnp <= 0.001  ~ "0.1%",
    #rnp <= 0.0001   ~ "0.01%",
    #rnp <= 0.00001  ~ "0.001%",
    TRUE ~ "else")
    ) ->
  collect.quant.glm.annot

collect.quant.glm.annot$rnp.thresh %>% table
####
#### -- join

rbind(mutate(collect.quant.glm.annot, subset = "glm"), 
      mutate(collect.quant.inv.annot, subset = "inv", rnp.thresh = "no.glm") ,fill = T) ->
  joint.qunt.data.for.plot

save(joint.qunt.data.for.plot, file = "fig2.fst.joint.qunt.data.for.plot.Rdata")

####
#### just cville
  city.select = c("Charlottesville")

  joint.qunt.data.for.plot %>%
  filter(pop %in% city.select) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    #color = type,
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Cville") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~as.factor(year_diff), scales =  "free_y") ->
  collect.quant.sum.cville

ggsave(collect.quant.sum.cville, file = "collect.quant.sum.cville.pdf", h = 5, w = 6)

###

joint.qunt.data.for.plot %>%
  filter(comp.set %in%  c("space.eu.e", "space.eu.w", "space.NoA.E") ) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    #color = type,
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Space sets in general") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~pop, scales =  "free_y") ->
  collect.quant.sum.space

ggsave(collect.quant.sum.space, file = "collect.quant.sum.space.pdf", h = 5, w = 6)

###
city.select.all = c("Akaa",  "Broggingen", "Cross Plains", "Linvilla", "Munich","Odesa" , "Yesiloz" )

joint.qunt.data.for.plot %>%
  filter(pop %in% city.select.all) %>%
  group_by(type, year_diff, pop, subset, rnp.thresh) %>%
  summarize(
    mean.p = mean(perc.above.control),
    sd.p = sd(perc.above.control),
    median.p = quantile(perc.above.control, 0.5),
    q25.p = quantile(perc.above.control, 0.25),
    q75.p = quantile(perc.above.control, 0.75)
  ) %>%
  ggplot(aes(
    x=type,
    y=median.p,
    ymin = q25.p,
    ymax= q75.p,
    shape = as.factor(year_diff),
    color =rnp.thresh
  )) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  geom_errorbar(width = 0.3, position = position_dodge(width = 0.5)) +
  geom_point( size = 2.1, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE) + 
  ggtitle("Cville") +
  ylab("Median of Percentiles (FST inv > FST controls) + IQRs") +
  xlab("Correlation to In(2L)t in DGRP") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  facet_grid(subset~pop, scales =  "free_y") ->
  collect.quant.sum.allcities

ggsave(collect.quant.sum.allcities, file = "collect.quant.sum.allcities.pdf", h = 9, w = 9)
