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

args <- commandArgs(trailingOnly = TRUE)

######
###### user defined parameters
##### objects needed for the analyses
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

load("/home/nzx3cc/snp_dt_25percMissing.Rdata")

setkey(snp.dt, id)

use <- apply(snp.dt[,c("VA_ch"), with=F],
             1, any)
snp.dt <- snp.dt[use]
setkey(snp.dt, id)

snp.dt <- snp.dt[N==0 & cm_mb>0 & !is.na(cm_mb) & chr%in%chrom]

#snps <- data.table(id = seqGetData(genofile, "variant.id"), miss = seqMissing(genofile))
#snps[id%in%snp.dt$id] %>% filter(miss <0.05) -> interest
#snp.dt %>% filter(id%in%interest$id) -> snp.dt
# 
# meta.addr <- "/home/nzx3cc/dest_v2.samps_26April2023.csv"
# 
# samps <- vroom(meta.addr)

#####$
#####
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
######
#filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))
# 
# 
# snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
#                       pos=seqGetData(genofile, "position"),
#                       variant.id=seqGetData(genofile, "variant.id"),
#                       nAlleles=seqNumAllele(genofile),
#                       missing=seqMissing(genofile), alleleFreq=seqAlleleFreq(genofile))


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



# snps.dt <- snps.dt[nAlleles==2]
# snps.dt <- snps.dt[alleleFreq > 0.05]
# snps.dt <- snps.dt[alleleFreq < 0.95]
# 
# snps.dt <- snps.dt[chr%in%c(chrom)][missing<.05]

if (thinningwindow > 0){
  startbp <- as.numeric(snp.dt$pos) + thinningwindow
  firstslice <- snp.dt[pos < startbp+thinningwindow,]
  pivotsnp <- firstslice[sample(.N, 1)]
  startbp <- pivotsnp[1, pos]
  thinnedsnps.dt <- data.table(pivotsnp)
  framenumber <- as.numeric(last(snp.dt$pos-first(snps.dt$pos))%/%thinningwindow)
  for (i in 1:framenumber) {
    thinnedframe <- snp.dt[pos > startbp & pos < (startbp + thinningwindow),]
    chosensnp <- data.table(thinnedframe[sample(.N,1)])
    startbp <- startbp + thinningwindow
    thinnedsnps.dt <- merge(thinnedsnps.dt, chosensnp, all = TRUE)
  }
thinnedsnps.dt -> snp.dt
  } 
#### ==
seqSetFilter(genofile,
             snp.dt$id)

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
  .$sampleId ->
  selected.samps
selected.samps <- unique(selected.samps)

ad.cville = ad[selected.samps, ]
dp.cville = dp[selected.samps, ]

####
####
<<<<<<< HEAD

write.table(snp.dt[,c(1,2,3,7)], paste0(outprefix, "_positiontable.txt"), quote=F)
=======
write.table(snp.dt[,c(1,2,3,7)], paste0(outprefix, "_positiontable.txt"))
>>>>>>> d32a99adb06c9f58a952da34379c36c7e094e371

pool <- new("pooldata",
            npools=dim(ad.cville)[1], #### Rows = Number of pools
            nsnp=dim(ad.cville)[2], ### Columns = Number of SNPs
            refallele.readcount=t(ad.cville),
            readcoverage=t(dp.cville),
            poolsizes=samps.cville$nFlies * 2,
            poolnames = rownames(ad.cville))
pooldata2genobaypass(pool, writing.dir = getwd(), prefix = outprefix)
