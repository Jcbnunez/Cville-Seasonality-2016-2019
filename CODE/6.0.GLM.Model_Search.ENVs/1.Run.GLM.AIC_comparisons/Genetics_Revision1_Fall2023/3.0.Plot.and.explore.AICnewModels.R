### Examine new models for Genetics revision #1

library(tidyverse)
library(data.table)
library(magrittr)

#### object 
newmods <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Revision_Best_Models/bestAIC.v2.Rdata"

  load(newmods)

###
 o2.ag <- o2.ag[!cluster%in%c("2.North_America_I95", "2.North_America_Mid", "2.North_America_W")]

### remake modRank to order by old-2L-inv

  focal_order <- o2.ag[perm_strategy=="old"][chr=="2L"][inv=="Inside Inversion"][cluster=="5.Cville"][,c("var", "mod", "modRank"), with=F]
  setnames(focal_order, "modRank", "modRank_focal")
  setkey(focal_order, var, mod)
  setkey(o2.ag, var, mod)
  o2.ag <- merge(o2.ag, focal_order)

### plot
  o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]
  o2.ag[var=="null", label:="null"]
  o2.ag[var=="pop_year", label:="pop_year"]
  o2.ag[is.na(label), label:="environment"]

  o2.ag[,perm_strategy:=factor(perm_strategy, levels=c("old", "new"))]
  o2.ag[,cluster:=factor(cluster, levels=c("5.Cville", "1.Europe_W", "3.Europe_E", "2.North_America_E"))]

  o2.ag[is.na(stage), stage:="stage0"]

  o2.ag.use <- o2.ag[stage=="stage0" | stage=="stage1" | (stage=="stage2" & !var%in%c("null", "pop_year"))]

####
  o2.ag.use[sig==T][perm_strategy=="new"] %>%
  group_by(chr, inv, cluster) %>%
  arrange(-rr) %>%
  slice_head(n = 10) ->
  best_top10_models

save(best_top10_models, file = "top10models.Rdata")
 
#### plot all
  p1 <- ggplot(data=o2.ag.use[sig==T][perm_strategy=="new"], aes(x=modRank_focal, y=log2(rr), 
  color=label, shape=label)) +
  facet_grid(cluster~chr+inv) +
  geom_point(data=o2.ag.use[sig==F][perm_strategy=="new"], aes(x=modRank_focal, y=log2(rr)), color="grey", alpha=.5) +
  geom_point() +
  geom_segment(data=o2.ag.use[var=="temp.max" & mod==2][stage=="stage2"],
              aes(x=modRank_focal+5, xend=modRank_focal, yend=log2(rr), y=-5), arrow = arrow(length = unit(0.1, "cm")), color="black") +
              theme(legend.pos = "bottom")


  ggsave(p1, file="./new_old.pdf", 
  width=12, height=7)

#### plot Virginia only

o2.ag.use[perm_strategy=="new"][cluster == "5.Cville"] -> dat.cville

dat.cville %>% 
filter(inv == "Inside Inversion",
chr == "2L"
) %>% arrange((rr)) %>%
mutate(in_rank = 1:nrow(.)) %>%
dplyr::select(var, mod, in_rank) ->
dat.cville.rank

dat.cville %>%
left_join(dat.cville.rank) ->
dat.cville.m

  cville.p <- ggplot(data=dat.cville.m[sig==T],
   aes(x=in_rank, y=log2(rr), 
  color=label, shape=label)) +
  facet_grid(inv~chr) + 
geom_hline(yintercept = 0) +  geom_point(data=dat.cville.m[sig==F], aes(x=in_rank, y=log2(rr)), color="grey", alpha=.5, size = 2.3) +
  geom_point(size = 2.3) +
  			  theme_bw() +
              theme(legend.pos = "bottom")

### ---> sig:=(prop.rr-2*prop.sd)>1

  ggsave(cville.p, file="./cville.p.pdf", 
  width=7, height=3)



