#SNP-wise analysis
# Verison 2 - August 19, 2022
# 
#
#

rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(gtools)
library(poolfstat)
library(viridis)

selected_pops = c("DE_Bro")
model = "humidity.ave;8;1.Europe_W"
popset="Brog"
####
####
####
####
####
####

### Load haplotags
load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplotag_snps_AFS_pol.Rdata")
haplotag_snps_AFS_pol %>% head
haplotag_snps_AFS_pol$win %>% table
#####
#####
#####
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
head(snp.dt)
#####
#####

#### import genomic data
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

samps %>%
  dplyr::select(samp1 = sampleId, locality_1 = locality,  year_1 = year, Month_1 = month, season_1 = season ) ->
  samp1

samps %>%
  dplyr::select(samp2 = sampleId, locality_2 = locality,  year_2 = year,  Month_2 = month,  season_2 = season ) ->
  samp2


### get subsample of data to work on
seqResetFilter(genofile)
#seqSetFilter(genofile, sample.id=samps.cville$sampleId)
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snp.dt$id)
snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]


#######
#######
#######
#######
#######

base <- "/project/berglandlab/alan/environmental_ombibus_global"
#########
file <- paste(base, model, 
                    paste(model, ".glmRNP.Rdata", sep = ""), 
                    sep = "/" )
print(file)
glm.out <- get(load(file))

glm.out %>%
  filter(chr == "2L") %>%
  filter(perm == 0) %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) -> 
  real.glm

#####
haplotag_snps_AFS_pol %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  .$SNP_id  %>% unique() -> snp.tags

haplotag_snps_AFS_pol %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(win %in% c("left", "right")) %>%
  group_by(win) %>%
  slice_head %>%
  .$SNP_id  %>% unique() -> inv.tags
  
tsp.tag = "2L_5192177_SNP"

#extract the haplotags
real.glm %<>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(SNP_id %in% c(snp.tags, inv.tags, tsp.tag) ) %>%
  mutate(win = case_when(
    pos > 2800000 & pos < 3200000 ~ "win_3.1",
    pos > 4470000 & pos < 4870000 ~ "win_4.7",
    pos > 4920000 & pos < 5320000 ~ "win_5.1",
    pos > 6000000 & pos < 6400000 ~ "win_6.1",
    pos > 6600000 & pos < 7000000 ~ "win_6.8",
    pos > 9400000 & pos < 9800000 ~ "win_9.6"
  )) 

#real.glm %>% filter(win %in% c("win_5.1") ) %>% arrange(rnp)

real.glm %>%
  filter(!is.na(win)) %>%
  group_by(win) %>%
  slice_min(rnp, with_ties = F) ->
  outliers.glm

rbind(filter(real.glm, SNP_id %in% c(inv.tags, tsp.tag)) , outliers.glm) ->
  outliers.glm

outliers.glm %<>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")])

####
#### --> select controls
glm.out %>%
  filter(chr != "2L") %>%
  filter(invName != "none") %>%
  filter(perm == 0) %>%
  filter(rnp > 0.25) %>%
  left_join(snp.dt[,c("chr", "pos", "cm_mb", "af")]) ->
  pool_of_controls

####
samps %>%
  filter(locality %in% selected_pops) ->
  samps.tmp

### get subsample of data to work on
seqResetFilter(genofile)

seqSetFilter(genofile, 
             sample.id=samps.tmp$sampleId, 
             variant.id=snp.dt$id)

dp.tmp <- seqGetData(genofile, "annotation/format/DP")

dp.tmp %>% t %>% as.data.frame() %>% 
  rowMeans(na.rm = T) -> per.snp.cov

data.frame(snp.dt[,c("chr", "pos")], cov.per.snp= per.snp.cov)-> cov.dt

####
in.df = left_join(outliers.glm, cov.dt)
cont.df = left_join(pool_of_controls, cov.dt)
chosen_controls = foreach(i = 1:dim(in.df)[1],
                          .combine = "rbind",
                          .errorhandling ="remove")%do%{
                            
  message(i)
    in.df[i,] ->
    tmp.achor
  
    cont.df %>%
    filter(cm_mb > tmp.achor$cm_mb-0.20 & cm_mb < tmp.achor$cm_mb+0.20,
           af > tmp.achor$af-0.030 & af < tmp.achor$af+0.030,
           cov.per.snp > tmp.achor$cov.per.snp-30 & cov.per.snp < tmp.achor$cov.per.snp+30
           ) %>%
    slice_sample(n = 200) %>%
    mutate(control_snp = paste(chr, pos, sep = "_" )) %>%
    mutate(matched_to = tmp.achor$SNP_id, type = "control")->
    select.control

  rbind(mutate(tmp.achor, matched_to = tmp.achor$SNP_id, type = "glm", SNP_id = paste(chr, pos, sep = "_") ), select.control, fill = T)
  
}

chosen_controls$matched_to %>% table

chosen_controls %<>% mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))
####
#### -----> Calculate FST
####
###
###
#samps %>% filter(continent == "NorthAmerica") %>% group_by(locality) %>% slice_head %>% 
#  select(city, locality) %>% as.data.frame()



######
#Generate outfile object

samps$Date = as.Date(samps.cville$Date, 
                            format='%Y/%m/%d')  

L = dim(samps.tmp)[1]

comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

print("Create combination vector")

comp_vector %<>%
  as.data.frame() %>%
  mutate(day_diff = NA, y_diff = NA)

##calculate day differences
print("Loop to esrtimate time difference")

for(i in 1:dim(comp_vector)[1]) {
  
  date1=samps.tmp$Date[comp_vector[i,1]]
  date2=samps.tmp$Date[comp_vector[i,2]]
  
  year1=samps.tmp$year[comp_vector[i,1]]
  year2=samps.tmp$year[comp_vector[i,2]]
  
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  comp_vector$y_diff[i] = abs(as.numeric(year1-year2))
  
  comp_vector$pop1[i] = samps.tmp$city[comp_vector[i,1]]
  comp_vector$pop2[i] = samps.tmp$city[comp_vector[i,2]]
  
  comp_vector$year1[i] = samps.tmp$year[comp_vector[i,1]]
  comp_vector$year2[i] = samps.tmp$year[comp_vector[i,2]]
  
  comp_vector$samp1[i] = samps.tmp$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = samps.tmp$sampleId[comp_vector[i,2]]
  
}

comp_vector %<>%
  mutate(city_test =
           ifelse(pop1 == pop2, "yes", "no"))

## only use observations from same population
print("Create comp vector")

comp_vector %>% 
  filter(y_diff %in% 0:1) %>%
  mutate(continent = 
           ifelse(.$pop1 %in% c("Charlottesville"),
                  "US", "EU" ))-> 
  comp_vector_samepops

if(popset=="cville"){
  comp_vector_samepops %<>%
    filter(year1 %in% 2016:2018 & year2 %in% 2015:2018) %>%
    filter(city_test == "yes")
}

print("End of part 3")
#######
#######

seqSetFilter(genofile, 
             sample.id=samps.tmp$sampleId, 
             variant.id=chosen_controls$variant.id)


### get allele frequency data
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), "SNP" ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , "SNP" , sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

colnames(dp) %in% sort(chosen_controls$SNP_id) %>% table
  
#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector_samepops)[1]),
  samp2 = rep(NA, dim(comp_vector_samepops)[1]),
  FST = rep(NA, dim(comp_vector_samepops)[1]),
  SNP.set = rep(NA, dim(comp_vector_samepops)[1])
)
snp.level.fst <- list()

for(i in 1:dim(comp_vector_samepops)[1]){
  
  print(i/dim(comp_vector_samepops)[1] * 100)
  
  samps_to_compare = c(comp_vector_samepops$samp1[i], comp_vector_samepops$samp2[i])
  
  pool_sizes = c(samps.tmp$nFlies[which(samps.tmp$sampleId == comp_vector_samepops$samp1[i])],
                 samps.tmp$nFlies[which(samps.tmp$sampleId == comp_vector_samepops$samp2[i])])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)
  
  
  fst.out <- computeFST(pool, method = "Anova")
  
  ### population level FST
  outfile$samp1[i] = comp_vector_samepops$samp1[i]
  outfile$samp2[i] = comp_vector_samepops$samp2[i]
  outfile$FST[i] = fst.out$FST
  outfile$SNP.set[i] = "glm.vs.control"
  
  #####################
  ### snp.level.fst
    fst.out$snp.FST %>% 
      data.frame(snp.FST = .) %>% 
      mutate(snp = rownames(.), 
             set = "glm.vs.control",
             samp1 = comp_vector_samepops$samp1[i],
             samp2 = comp_vector_samepops$samp2[i]) -> tmp.fst.snp.wise
  
    snp.level.fst[[i]] = tmp.fst.snp.wise
    
}##  close i   

### save snp.wise.object
snp.wise.tmp.f.df = do.call(rbind.data.frame, snp.level.fst)

snp.wise.tmp.f.df %>%
  left_join(samp1) %>%
  left_join(samp2) %>% 
  mutate(dy = abs(year_1-year_2)) ->
  snp.wise.tmp.f.df.raw

save(snp.wise.tmp.f.df.raw, 
     file = paste(popset,"RAW.snp.wise.df.Rdata", sep = "." ))

snp.wise.tmp.f.df.raw %>%
  #filter(year_1 %in% 2016:2018 & year_2 %in% 2016:2018) %>%
  group_by(SNP_id = snp, dy) %>%
  summarize(mean.fst = mean(snp.FST, na.rm = T),
            median.fst = median(snp.FST, na.rm = T)) %>%
  left_join(chosen_controls) -> snp.wise.tmp.f.df.ag

save(snp.wise.tmp.f.df.ag, 
     file = paste(popset,"AG.snp.wise.df.Rdata", sep = "." ))

####
dat.in.qtile = filter(snp.wise.tmp.f.df.ag, type == "glm", dy %in% 0:1)
controls.qtile = filter(snp.wise.tmp.f.df.ag, type == "control", dy %in% 0:1)

quantiles = foreach(i=1:dim(dat.in.qtile)[1], .combine = "rbind")%do%{
  
  dy.set = dat.in.qtile$dy[i]
  SNP_id.set = dat.in.qtile$SNP_id[i]
  fst.glm = dat.in.qtile$median.fst[i]
  
  message(paste(i , dy.set, SNP_id.set, fst.glm,  sep = " | "))
  
  controls.qtile %>%
    filter(dy == dy.set & matched_to == SNP_id.set ) ->
    fst.dist.tmp
  
  quant.glm = sum(sort(fst.dist.tmp$median.fst) > fst.glm)/length(fst.dist.tmp$median.fst)
  
  data.frame(matched_to = SNP_id.set, 
             dy =  dy.set, 
             quntile =  quant.glm,
             pop = popset)
  
}

save(quantiles, 
     file = paste(popset,"quantiles.Rdata", sep = "." ))

quantiles %>%
  group_by(dy, pop) %>%
  summarize(mean(quntile))

####
####


ggplot() +
  geom_boxplot(data = filter(snp.wise.tmp.f.df.ag, type == "control", dy %in% 0:1
  ),
  aes(x=as.factor(dy),
      color = as.factor(dy),
      y=median.fst),
  outlier.shape = NA, width = 0.6) +
  geom_point(data = left_join(filter(snp.wise.tmp.f.df.ag, type == "glm", dy %in% 0:1), quantiles),
             aes(x=as.factor(dy),
                 y=median.fst),
             size = 5,
             color = "red") +
  geom_text(data = left_join(filter(snp.wise.tmp.f.df.ag, type == "glm", dy %in% 0:1), quantiles),
             aes(x=as.factor(dy),
                 y=median.fst,
                 label = round(quntile, 2 )),
             size = 2) +
  ylab("Median FST (SNP-wise)") +
  ggtitle("Brog") +
  #ylim(-0.02, 0.005) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 6),
        legend.pos = "bottom") +
  facet_grid(.~matched_to) ->
  box.dy.fst

ggsave(box.dy.fst, file = paste(popset, "box.dy.fst.pdf", sep = "."), w = 9, h = 4)
