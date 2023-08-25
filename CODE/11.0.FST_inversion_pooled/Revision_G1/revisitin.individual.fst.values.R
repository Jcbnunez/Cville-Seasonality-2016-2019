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

# "/sfs/weka/scratch/yey2sn/Genetics_Resubmission/19.inv.fst" 
rm(list=ls())

### This gets us individual SNP fsts
load("cville.V3.snp.wise.df.Rdata")
setDT(snp.wise.tmp.f.df)

#### info about glms
glm_info <- get(load("Cville.GLM.rnp5.snps.Rdata"))
glm_info %>%
dplyr::select(snp=SNP_id, rnp) ->
glm_info


#########

#### load
cont.info <- fread("Matched.controls.In2Lt.txt")
setDT(cont.info)

cont.info$SNP_id %>% unique -> controls
cont.info$matched_to %>% unique -> glms

snp.wise.tmp.f.df %>%
filter(snp %in% controls)


###
snp.wise.tmp.f.df %>% 
separate(snp, remove = F, into = c("chr","pos","fetaure"), sep = "_") %>%
mutate(snp_type = 
case_when(
chr == "2L" ~ "glms",
TRUE ~ "controls"
)) %>%
  mutate(WinSpec = case_when(
    #chr == "2L" & pos > 2000000 & pos < 2400000 ~ "L",
    #chr == "2L" & pos > 12900000 & pos < 13300000 ~ "R",
    #chr == "2L" & pos > 2800000 & pos < 3200000 ~ "W3.1",
    #chr == "2L" & pos > 4470000 & pos < 4870000 ~ "W4.6",
    #chr == "2L" & pos > 4920000 & pos < 5320000 ~ "W5.2",
    #chr == "2L" & pos > 6000000 & pos < 6400000 ~ "W6.1",
    #chr == "2L" & pos > 6600000 & pos < 7000000 ~ "W6.8",
    #chr == "2L" & pos > 9400000 & pos < 9800000 ~ "W9.6",
    chr == "2L" ~ "glm",
    TRUE ~ "NW")) -> snp.wise.tmp.f.df.annot

snp.wise.tmp.f.df.annot %>%
group_by(WinSpec) %>%
summarize(mfst = mean(snp.FST, na.rm = T))


### Add some info for temperature and time
load("./DatFor.Haplotypes.trajectory.time.weather.Rdata")
Cville_haplotags_for_viz %>%
dplyr::select(sampleId, temp.max) %>%
group_by(sampleId) %>%
slice_head() ->
weather.dat

weather.dat1=weather.dat
names(weather.dat1) = c("samp1","weth1")
weather.dat2=weather.dat
names(weather.dat2) = c("samp2","weth2")

### more stuff
master <- fread("master.comp.list.fst.txt")
master %>%
filter(year_diff == 0) %>%
dplyr::select(day_diff, year_diff, samp1, samp2) ->
master.sub


##
log10_ceiling <- function(x) {
    10^(ceiling(log10(x)))
}
##

####
snp.wise.tmp.f.df.annot %>%
left_join(weather.dat1) %>%
left_join(weather.dat2) %>%
mutate(deltaT = abs(weth1-weth2))  %>% 
left_join(master.sub) %>%
full_join(glm_info) %>% 
mutate(rnp.ceil = log10_ceiling(rnp)) ->
snp.wise.tmp.f.df.annot.weather


o =
foreach(i = c("W3.1" , "W4.6", "W5.2" ,"W6.8", "W9.6"),
.combine = "rbind")%do%{
summary(lm(snp.FST ~ day_diff, 
data = filter(snp.wise.tmp.f.df.annot.weather,
WinSpec == i,
chr == "2L"))) -> out
data.frame(i, out$coefficients[2,1], out$coefficients[2,4])

}

summary(lm(snp.FST ~ day_diff, 
data = filter(snp.wise.tmp.f.df.annot.weather,
rnp.ceil == 1e-04,
chr == "2L"))) 
summary(lm(snp.FST ~ day_diff, 
data = filter(snp.wise.tmp.f.df.annot.weather,
rnp.ceil == 1e-03,
chr == "2L"))) 


oC =
foreach(i = c("NW"),
.combine = "rbind")%do%{
summary(lm(snp.FST ~ day_diff, 
data = filter(snp.wise.tmp.f.df.annot.weather,
WinSpec == i,
chr != "2L", year_diff == 0))) -> out
data.frame(i, out$coefficients[2,1], out$coefficients[2,4])

}


## Temperature
o2 =
foreach(i = c("W3.1" , "W4.6", "W5.2" ,"W6.8", "W9.6"),
.combine = "rbind")%do%{

summary(lm(snp.FST ~ deltaT, 
data = filter(snp.wise.tmp.f.df.annot.weather,
WinSpec == i,
chr == "2L"))) -> out

data.frame(i, out$coefficients[2,1], out$coefficients[2,4])

}
#not in window -- controls
summary(lm(snp.FST ~ deltaT, 
data = filter(snp.wise.tmp.f.df.annot.weather,
WinSpec == "NW",
chr != "2L"))) 





cor.test(~ snp.FST + deltaT, 
data = filter(snp.wise.tmp.f.df.annot.weather,
WinSpec == "W3,1",
chr == "2L"))

###
snp.wise.tmp.f.df.annot.weather %>%
filter(!is.na(day_diff)) %>%
ggplot(aes(
x=day_diff,
y=snp.FST,
color = snp_type
)) +
#geom_point() +
geom_smooth(method = "lm") +
facet_grid(~paste(chr,WinSpec)) ->
plot.dt.fst
ggsave(plot.dt.fst, file = "plot.dt.fst.pdf", w = 12, h = 3.5)


### measure
snp.wise.tmp.f.df.annot %>%
group_by(snp_type) %>%
summarize(m.fst = mean(snp.FST, na.rm = T))

##
log10_ceiling <- function(x) {
    10^(ceiling(log10(x)))
}
##

snp.wise.tmp.f.df.annot %>%
full_join(glm_info) %>% 
mutate(rnp.ceil = log10_ceiling(rnp)) %>% 
group_by(snp_type, rnp.ceil ) %>%
summarize(m.fst = mean(snp.FST, na.rm = T))

###
#  snp_type rnp.ceil      m.fst
#  <chr>       <dbl>      <dbl>
#1 controls NA        -0.000426
#2 glms      0.0001    0.0153  
#3 glms      0.001     0.00948 

