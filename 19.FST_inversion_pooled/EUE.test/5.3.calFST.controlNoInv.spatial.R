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
library(geosphere)

#### load snp metadata
load("./EUE.glm.and.matchedControls.Rdata")
set = "macthed.controls.noInv"
cluster =  "3.Europe_E"
anlysis = "geo"

#glm.snps
#macthed.controls.Inv
#macthed.controls.noInv
selected_Set = macthed.controls.noInv


####

#### import genomic data
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]


###
### import demo clusters
file.clusters = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/DEST_Sample_clusters.txt"
clust.asg = fread(file.clusters)

clust.asg %>%
  filter(Continental_clusters == cluster ) ->
  clust.dat

samps %>%
  filter(set == "DrosEU") %>%
  filter(sampleId %in% clust.dat$sampleId ) ->
  samps.cville

samps.cville %<>%
  filter(year == 2015) %>%
  filter(season == "fall")


### get subsample of data to work on
seqResetFilter(genofile)

seqSetFilter(genofile, 
             sample.id=samps.cville$sampleId, 
             variant.id=selected_Set$id)

######
#Generate outfile object

samps.cville$Date = as.Date(samps.cville$Date, 
                             format='%Y/%m/%d')  

L = dim(samps.cville)[1]

comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

print("Create combination vector")

comp_vector %<>%
  as.data.frame() %>%
  mutate(day_diff = NA)

##calculate day differences
print("Loop to esrtimate time difference")

for(i in 1:dim(comp_vector)[1]) {
  
  date1=samps.cville$Date[comp_vector[i,1]]
  date2=samps.cville$Date[comp_vector[i,2]]
  
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  
  comp_vector$pop1[i] = samps.cville$city[comp_vector[i,1]]
  comp_vector$pop2[i] = samps.cville$city[comp_vector[i,2]]
  
  comp_vector$samp1[i] = samps.cville$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = samps.cville$sampleId[comp_vector[i,2]]
  
}

comp_vector %<>%
  mutate(city_test =
           ifelse(pop1 == pop2, "yes", "no"))

## only use observations from same population
print("Create comp vector")

comp_vector %>% 
  .[which(.$city_test == "no"),] %>%
  #.[which(.$day_diff < 360),] %>%
  mutate(continent = 
           ifelse(.$pop1 %in% c("Charlottesville"),
                  "US", "EU" ))-> 
  comp_vector_samepops

## add lat long
samps.cville %>%
  dplyr::select(samp1 = sampleId, lat1=lat, long1=long) -> met1

samps.cville %>%
  dplyr::select(samp2 = sampleId, lat2=lat, long2=long) -> met2

comp_vector_samepops %<>%
  left_join(met1) %>%
  left_join(met2)

comp_vector_samepops %<>%
  group_by(samp1, samp2) %>%
  mutate(geo.dist = distVincentyEllipsoid(c(long1,lat1),
                                          c(long2,lat2))) 

print("End of part 3")
#######
#######

### get allele frequency data
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")


#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector_samepops)[1]),
  samp2 = rep(NA, dim(comp_vector_samepops)[1]),
  FST = rep(NA, dim(comp_vector_samepops)[1]),
  SNP.set = rep(NA, dim(comp_vector_samepops)[1])
)


for(i in 1:dim(comp_vector_samepops)[1]){
  
  print(i/dim(comp_vector_samepops)[1] * 100)
  
  samps_to_compare = c(comp_vector_samepops$samp1[i], comp_vector_samepops$samp2[i])
  
  pool_sizes = c(samps.cville$nFlies[which(samps.cville$sampleId == comp_vector_samepops$samp1[i])],
                 samps.cville$nFlies[which(samps.cville$sampleId == comp_vector_samepops$samp2[i])])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)
  
  
  fst.out <- computeFST(pool, method = "Anova")
  
  outfile$samp1[i] = comp_vector_samepops$samp1[i]
  outfile$samp2[i] = comp_vector_samepops$samp2[i]
  outfile$FST[i] = fst.out$FST
  outfile$SNP.set[i] = set
  
}##  close i   

#### end fst calc, next step below
left_join(comp_vector_samepops, outfile) -> Out_comp_vector_samepops

samps %<>%
  mutate(month_col = month(as.Date(collectionDate, 
                                   format = c("%m/%d/%Y"))))

samps %>%
  dplyr::select(sampleId, month_col, year) -> samp_1_meta
names(samp_1_meta) = c("samp1", "month1", "year1")

samps %>%
  dplyr::select(sampleId, month_col, year) -> samp_2_meta
names(samp_2_meta) =  c("samp2", "month2", "year2")

Out_comp_vector_samepops %<>%
  left_join(samp_1_meta) %>% 
  left_join(samp_2_meta) 

print("Save object")
save(Out_comp_vector_samepops, 
     file = paste(set,anlysis,"EUE.Year_to_year_object.Rdata", sep = "." ))
