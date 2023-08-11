rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(magrittr)
library(reshape2)
library(ape)
library(patchwork)


### create array of files
fold <- "/gpfs2/scratch/jcnunez/genetics_resub/2.spatial_fst/spatial.fst"

filslod =
system(paste("ls ",fold, "/*", sep = ""), intern = T)

### load

infst <-
foreach(i = filslod, 
.combine = "rbind")%do%{

message(i)
tmp <- get(load(i))

}

setDT(infst)
### some filtering
infst %>%
separate(samp1,
into = c("country1","city1","year1","id1"), sep = "_"
, remove = F) %>%
separate(samp2,
into = c("country2","city2","year2","id2"), sep = "_", 
remove = F) %>% 
mutate(loc_comp = paste(
paste(country1,city1, sep = "_"),
paste(country2,city2, sep = "_"), 
sep = "|")) -> infst.mod

### exploration
infst.mod %>%
dcast(loc_comp~year, value.var = "FST") %>%
.[complete.cases(.)] ->
fst.matrix.spatial

#### summarize
infst.mod %>%
filter(loc_comp %in% fst.matrix.spatial$loc_comp) ->
symetrical_pairs
unique(c(symetrical_pairs$samp1, symetrical_pairs$samp2))

#### heat plot
fst.matrix.spatial %>%
ggplot(aes(
x=`2014`,
y=`2015`)) +
geom_point() +
geom_abline(slope = 1) ->
cor1415

fst.matrix.spatial %>%
ggplot(aes(
x=`2015`,
y=`2016`)) +
geom_point() +
geom_abline(slope = 1) ->
cor1516

fst.matrix.spatial %>%
ggplot(aes(
x=`2014`,
y=`2016`)) +
geom_point() +
geom_abline(slope = 1) ->
cor1416

ggsave(cor1415+cor1516+cor1416, file = "corr.fst.pdf", h = 3, w = 9)

#### correlations
cor.test(~`2014`+`2015`, data = fst.matrix.spatial)
cor.test(~`2015`+`2016`, data = fst.matrix.spatial)
cor.test(~`2014`+`2016`, data = fst.matrix.spatial)

#### summarize
symetrical_pairs %>%
ggplot(
aes(
x=hav_d,
y=FST,
color = as.factor(year)
)) +
geom_point() +
geom_smooth(method = "lm") +
facet_wrap(~cont1) ->
fst.havd

ggsave(fst.havd, file = "fst.havd.pdf", h = 3, w = 7)

summary(lm(FST ~ hav_d +year, data = symetrical_pairs))
