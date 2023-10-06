library(tidyverse)
library(magrittr)
library(vroom)
library(lubridate)
library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)

library(foreach)
library(doMC)
library(doParallel)
registerDoMC(10)



load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;5;5.Cville/temp.max;5;5.Cville.glmRNP.Rdata")
glm.out <- glm.out[perm==0]
glm.out -> snp.dt
args <- commandArgs(trailingOnly = TRUE)
print(args)
chrom <- as.character(args[2])
if (chrom == "all"){
  chrom <- c("2L", "2R", "3L", "3R")
} else if (chrom =="no2L"){
  chrom <- c("2R", "3L", "3R")
} else if (chrom == "no2R"){
  chrom <- c("2L", "3L", "3R")
} else if (chrom == "no3L"){
  chrom <- c("2L", "2R", "3R")
} else if (chrom == "no3R"){
  chrom <- c("2L", "2R", "3L")
}
outprefix <- args[3]
print(typeof(outprefix))
print(outprefix)
thinningwindow <- as.numeric(args[4])

genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

samps[locality=="UA_od", locality:="UA_Ode"]
weather.ave <- read.table("/scratch/nzx3cc/nzx3cc/env_analysis/weather.ave")
samps <- merge(samps, weather.ave, by="sampleId")

### define five population sets: West-Coast NA; East-Coast NA; Western Europe; Eastern Europe; Cville
### also  clean up based on other flags from Kapun et al
dest <- fread("https://raw.githubusercontent.com/DEST-bio/DEST_freeze1/main/populationInfo/samps_10Nov2020.csv")
samps <- merge(samps, dest[status=="Keep"][propSimNorm<.05][,"sampleId"], by="sampleId", all.x=T)
cluster <- fread("https://raw.githubusercontent.com/DEST-bio/DEST_freeze1/main/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")
samps <- merge(samps, cluster[,c("sampleId", "Continental_clusters"), with=F], by="sampleId", all.x=T)
samps[,cluster:=Continental_clusters]

samps[long< -110, cluster:="2.North_America_W"]
samps[long> -110 & Continental_clusters=="2.North_America", cluster:="2.North_America_E"]
samps[locality=="VA_ch", cluster:="5.Cville"]

samps[locality%in%c("IL_ur", "KA_to", "MI_bh", "ON_su", "WI_cp"), cluster:="2.North_America_Mid"]

setkey(samps, sampleId)
table(samps[!duplicated(samps)]$cluster)



if (thinningwindow > 0){
  startbp <- as.numeric(snp.dt$pos) + thinningwindow
  firstslice <- snp.dt[pos < startbp+thinningwindow,]
  pivotsnp <- firstslice[sample(.N, 1)]
  startbp <- pivotsnp[1, pos]
  thinnedsnps.dt <- data.table(pivotsnp)
  framenumber <- as.numeric(last(snp.dt$pos-first(snp.dt$pos))%/%thinningwindow)
  for (i in 1:framenumber) {
    thinnedframe <- snp.dt[pos > startbp & pos < (startbp + thinningwindow),]
    chosensnp <- data.table(thinnedframe[sample(.N,1)])
    startbp <- startbp + thinningwindow
    thinnedsnps.dt <- rbind(thinnedsnps.dt, chosensnp)
  }
  thinnedsnps.dt -> snp.dt
} 
#### ==
seqSetFilter(genofile,
             snp.dt$variant.id)

#### ==>

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
#### Filter arbitrarily
samps %>%
  filter(set == "CvilleSet") ->
  samps.cville

samps.cville %>%
  .$sampleId -> selected.samps
selected.samps <- unique(selected.samps)
ad.cville = ad[selected.samps, ]
dp.cville = dp[selected.samps, ]
samps.cville$nFlies -> nflies
####
####
write.table(snp.dt[,c("pos","chr","variant.id","invName")], paste0(outprefix, "_positiontable.txt"), quote=F)
pool <- new("pooldata",
            npools=dim(ad.cville)[1], #### Rows = Number of pools
            nsnp=dim(ad.cville)[2], ### Columns = Number of SNPs
            refallele.readcount=t(ad.cville),
            readcoverage=t(dp.cville),
            poolsizes=nflies * 2,
            poolnames = rownames(ad.cville))
pooldata2genobaypass(pool, writing.dir = getwd(), prefix = outprefix)

