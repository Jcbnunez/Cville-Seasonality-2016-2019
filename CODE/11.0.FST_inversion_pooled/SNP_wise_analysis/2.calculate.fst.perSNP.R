##### --- Calculate FST 
##### 

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

### ---> user defined parameter
###array.job=1
args = commandArgs(trailingOnly=TRUE)
array.job=as.numeric(args[1])

#########
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

#####
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
#head(snp.dt)
snp.dt %>% dim %>% .[1] -> max.snp.num
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=snp.dt$id)
snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

win.bp <- 500000
step.bp <- 500000+1

setkey(snp.dt, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- snp.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
## --> 191 jobs

### create the snp dt parition

 snp.dt %>%
   filter(chr == wins$chr[array.job] ) %>%
   filter(pos >= wins$start[array.job] & pos <= wins$end[array.job]  ) ->
   snp.dt.tmp
 
#####
comparisons <- fread("./master.comp.list.fst.txt")
#head(comparisons)
#####

# over comparisons

snps.fst.dt = foreach(comp = 1:dim(comparisons)[1], 
        .combine = "rbind", 
        .errorhandling = "remove")%do%{
  message(comp)
  
  comparisons %>%
        .[comp,] ->
        comp.tmp
          
  samps.tmp = c(comp.tmp$samp1, comp.tmp$samp2)
  
    ### get subsample of data to work on
    seqResetFilter(genofile)
    #seqSetFilter(genofile, sample.id=samps.cville$sampleId)
    seqSetFilter(genofile, sample.id=samps.tmp, 
                 variant.id=snp.dt.tmp$id)
    
    print("Create ad and dp objects")
    
    ad <- seqGetData(genofile, "annotation/format/AD")
    ad <- ad$data
    dp <- seqGetData(genofile, "annotation/format/DP")
    
    print("Create dat object")
    #message(paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), "SNP" ,  sep="_"))
    
    #Add metadata ad
    colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), "SNP" ,  sep="_")
    rownames(ad) <- seqGetData(genofile, "sample.id")
    
    #Add metadata dp
    colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , "SNP" , sep="_")
    rownames(dp) <- seqGetData(genofile, "sample.id")
    
    samps_to_compare = c(comp.tmp$samp1, comp.tmp$samp2)
    
    pool_sizes = c(samps$nFlies[which(samps$sampleId == comp.tmp$samp1)],
                   samps$nFlies[which(samps$sampleId == comp.tmp$samp2)])
    
    
    ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
    rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
    
    pool <- new("pooldata",
                npools=dim(ad.matrix)[1], #### Rows = Number of pools
                nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
                refallele.readcount=t(ad.matrix),
                readcoverage=t(rd.matrix),
                poolsizes=pool_sizes * 2)
    
    fst.out <- computeFST(pool, method = "Anova")
    
    data.frame(SNP_id = names(fst.out$snp.FST), SNPwise.FST = fst.out$snp.FST) -> snpfst
    data.frame(SNP_id = names(fst.out$snp.FST), SNPwise.Q1 = fst.out$snp.Q1) -> snpq1
    data.frame(SNP_id = names(fst.out$snp.FST), SNPwise.Q2 = fst.out$snp.Q2) -> snpq2
    
    
    data.frame(
      comp.tmp,
      snp.dt.tmp
    ) %>%
      mutate(
        window.FST = fst.out$FST,
        window.name = paste(wins$chr[array.job],wins$start[array.job], 
                            wins$end[array.job], sep = "_"),
        window.mid.pt = (wins$start[array.job]+wins$end[array.job])/2
        ) -> metadat
    
    tmp.fst.out =
    cbind(metadat, snpfst, SNPwise.Q1=snpq1[,-1], SNPwise.Q2=snpq2[,-1])

    return(tmp.fst.out)
}


 snps.fst.dt %>%
   .[complete.cases(.$SNPwise.FST),] %>%
   mutate(job.array.id = array.job) ->
   snps.fst.dt.flt

 ### save
 
out.folder = "/project/berglandlab/DEST_fst_out/out.files/"
naming.convention = paste("SNPFST",
                          array.job, 
                          wins$chr[array.job],wins$start[array.job], 
                          wins$end[array.job], sep = "_" )

save(snps.fst.dt.flt,
     file = paste(out.folder, naming.convention, ".Rdata", sep = "" ))
