library(tidyverse)
library(magrittr)
library(vroom)

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
  chrom <- c("2L", "2R", "3L", "3R")}
outprefix <- args[3]
print(typeof(outprefix))
print(outprefix)
thinningwindow <- as.numeric(args[4])

meta.addr <- "/home/nzx3cc/dest_v2.samps_26April2023.csv"

samps <- vroom(meta.addr)

#####$
#####
genofile.path <- "/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

######
#filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))


snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile), alleleFreq=seqAlleleFreq(genofile))

#snp.dt <- snp.dt[cm_mb>0 & !is.na(cm_mb)]

snps.dt <- snps.dt[nAlleles==2]
snps.dt <- snps.dt[alleleFreq > 0.05]
snps.dt <- snps.dt[alleleFreq < 0.95]

snps.dt <- snps.dt[chr%in%c(chrom)][missing<.05]

if (thinningwindow > 0){
startbp <- as.numeric(snps.dt[1, 2]) + thinningwindow
firstslice <- snps.dt[pos < startbp+thinningwindow,]
pivotsnp <- firstslice[sample(.N, 1)]
startbp <- pivotsnp[1, pos]
thinnedsnps.dt <- data.table(pivotsnp)
framenumber <- (as.numeric(last(snps.dt)[,2]-first(snps.dt)[,2]))%/%thinningwindow
for (i in 1:framenumber) {
thinnedframe <- snps.dt[pos > startbp & pos < (startbp + thinningwindow),]
chosensnp <- data.table(thinnedframe[sample(.N,1)])
startbp <- startbp + thinningwindow
thinnedsnps.dt <- merge(thinnedsnps.dt, chosensnp, all = TRUE)
}
} else {
  snps.dt -> thinnedsnps.dt
}
#### ==
seqSetFilter(genofile,
             thinnedsnps.dt$variant.id)

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
  filter(set == "cville") ->
  samps.cville

samps.cville %>%
  .$sampleId ->
  selected.samps

ad.cville = ad[selected.samps, ]
dp.cville = dp[selected.samps, ]

####
####

pool <- new("pooldata",
            npools=dim(ad.cville)[1], #### Rows = Number of pools
            nsnp=dim(ad.cville)[2], ### Columns = Number of SNPs
            refallele.readcount=t(ad.cville),
            readcoverage=t(dp.cville),
            poolsizes=samps.cville$nFlies * 2,
            poolnames = rownames(ad.cville))
pooldata2genobaypass(pool, writing.dir = getwd(), prefix = outprefix)
