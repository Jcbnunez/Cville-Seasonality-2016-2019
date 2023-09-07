
### library
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(foreach)

###
glm.object <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Revision_Best_Models/temp.max;2;5.Cville.v2.glmRNP.Rdata"

load(glm.object)

Inv.NoInv.bin =
foreach( p = 0:100, .combine = "rbind")%do%{
foreach( ch = c("2L","2R","3R","3L"), .combine = "rbind")%do%{
foreach( invName = c("Inv","noInv"), 
.combine = "rbind")%do%{

message(paste(p, ch, invName, sep = "*"))

if(invName == "Inv"){
glm.out %>%
filter(perm == p) %>%
filter(chr == ch) %>%
filter(invName != "none" ) %>%
.$p_lrt.x %>%
hist(breaks = 100) -> hist.obj 

data.frame(
hist.obj$mids,
hist.obj$counts,
perm = p,
chr = ch,
inv = "Inv"
) -> o
}

if(invName == "noInv"){
glm.out %>%
filter(perm == p) %>%
filter(chr == ch) %>%
filter(invName == "none" ) %>%
.$p_lrt.y %>%
hist(breaks = 100) -> hist.obj 

data.frame(
hist.obj$mids,
hist.obj$counts,
perm = p,
chr = ch,
inv = "NoInv"
) -> o
}

return(o)
}}}

save(Inv.NoInv.bin, file = "Inv.NoInv.bin.Rdata")

###

Inv.NoInv.bin %>%
ggplot(aes(
x=(hist.obj.mids),
y=hist.obj.counts,
group=as.factor(perm),
color = perm==0,
size = perm==0
)) +
geom_line(aes(alpha =perm==0 )) +
scale_alpha_manual(values = c(0.1, 1)) +
scale_size_manual(values = c(0.7, 1.3)) +
scale_color_manual(values = c("pink","black")) +
theme_bw() +
scale_x_log10() +
facet_grid(inv~chr)->
bin.plot

ggsave(bin.plot, file = "bin.pdf", w = 8, h = 3)

### data --- Alans code:
#  setwd("/project/berglandlab/alan/environmental_ombibus_global_permV2/glmBins")

#=#enrich <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Revision_Best_Models/Enrichment_objects/temp.max;2;5.Cville.v2.glm_enrich.Rdata"
#=#
#=#  load(enrich)
#=#
#=#### plot
#=#glm.enrich %>%
#=#ggplot(aes(
#=#color=perm==0,
#=#x=(thr),
#=#y=-log10(new)
#=#)) +
#=#geom_line(aes(group=as.factor(perm))) +
#=#facet_wrap(invName=="none"~chr, scales = "free") ->
#=#bins.enrich
#=#
#=#ggsave(bins.enrich, file = "bins.enrich.pdf")
#=#
#=#
#=##  glm.bins[,list(n=sum(N)), list(perm, perm_method)]
#=#  ggplot(data=glm.bins[perm<10], aes(x=-log10(p_lrt.low), y=N, group=perm, color=as.factor(perm!=0))) + geom_line() + facet_grid(perm_method~chr+inv)
