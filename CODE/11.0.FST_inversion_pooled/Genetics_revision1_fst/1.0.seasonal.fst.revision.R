rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(magrittr)
library(SeqArray)
library(lubridate)
library(gtools)
library(poolfstat)
library(geosphere)

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1]) # 1-469
### Shift is a parameter to add to k, to compensate for the array limit of VACC
#shift=as.numeric(args[2])

#print(kor)
#print(shift)

#k=kor+shift
print(k)


#load genofile
genofile <- seqOpen("../dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

#samps
load("../DEST_EC_metadata.Rdata")
samps = samps_EFFCOV %>% 
filter(MeanEC >= 28) %>% 
filter(set == "CvilleSet")

### Add temperature
#system("cp /netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/weatherAve.Rdata ../")
#load("../weatherAve.Rdata")
load("../DatFor.Haplotypes.trajectory.time.weather.Rdata")

Cville_haplotags_for_viz %>%
dplyr::select(sampleId, temp.max) %>%
group_by(sampleId) %>%
slice_head() ->
weather.dat

####
left_join(samps, weather.dat) ->
samps.weather

####
cville.differential =
foreach(y = 2016:2018, .combine = "rbind")%do%{

tmp.samps = samps.weather %>%
filter(year == y)

message(y)

tmp.samps$collectionDate = as.Date(tmp.samps$collectionDate, format='%m/%d/%Y')  

L = dim(tmp.samps)[1]

comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

comp.vector =
foreach(i = 1:dim(comp_vector)[1],
.combine = "rbind")%do%{

  message(i)

  samp1 = tmp.samps$sampleId[comp_vector[i,1]]
  samp2 = tmp.samps$sampleId[comp_vector[i,2]]
 
  pool1 = tmp.samps$nFlies[comp_vector[i,1]]
  pool2 = tmp.samps$nFlies[comp_vector[i,2]]
 
  Temp1 = tmp.samps$temp.max[comp_vector[i,1]]
  Temp2 = tmp.samps$temp.max[comp_vector[i,2]]

 
  data.frame(
  samp1,
  samp2,
  pool1,
  pool2,
  Temp1,
  Temp2,
  year=y
  ) %>%
   mutate(Tdelta = abs(Temp1-Temp2) )
  
} ## i

return(comp.vector)

}

#### create differentials
cville.differential %<>%
mutate(
delta_T_group =
case_when(
Tdelta < 5 ~ "dT<5",
Tdelta >= 5 & Tdelta <= 10 ~ "5<dT<10",
Tdelta > 10 ~ "dT>10",
))

####

## filter genoefile
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L")]
#,"2R","3L","3R"

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

##### PREPARE WINDOWING ALGORITHM
  win.bp <- 1e5
  step.bp <- 5e4
  
  setkey(snps.dt , "chr")
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L"
  					#,"2R", "3L", "3R"
  					),
                  .combine="rbind", 
                  .errorhandling="remove")%do%{
                    
                    tmp <- snps.dt %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  
  wins[,i:=1:dim(wins)[1]]
  
  dim(wins)
######


###### ----> filter by windows
  setkey(snps.dt, chr, pos)


  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind"
  )%do%{

    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.snp.dt <- snps.dt[J(data.table(chr=wins[win.i]$chr,                                 pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

chr=wins[win.i]$chr
start.p=wins[win.i]$start
end.p=wins[win.i]$end

seqResetFilter(genofile)

seqSetFilter(genofile, 
sample.id=samps$sampleId,
variant.id=win.snp.dt$variant.id)

## create count -- coverage objects
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

## add metadata
sampleids <- seqGetData(genofile, "sample.id")

colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")
####

guide_in = cville.differential

#o =
#foreach(
#k = 1:dim(guide_in)[1] ,
# .combine = "rbind",
# .errorhandling = "remove")%do%{

#print(k/dim(guide_in)[1] * 100)

#k = 1
##### ===>>> K is declared above

  samps_to_compare = c(guide_in$samp1[k], guide_in$samp2[k])
  pool_sizes = c(guide_in$pool1[k], guide_in$pool2[k])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  

  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)

  fst.out <- computeFST(pool, method = "Anova")
  
  guide_in[k,] %>%
  mutate(FST = fst.out$FST,
  chr=wins[win.i]$chr,
  start=wins[win.i]$start,
  end=wins[win.i]$end
  ) ->
  o.tmp
 
# return(o.tmp)

#}

#####

return(o.tmp)

}

####

root = "/gpfs2/scratch/jcnunez/genetics_resub/3.seasonal_fst/out_file/"

save(win.out,
file = paste(root, "seasFST." ,k, ".Rdata", sep = ""))

