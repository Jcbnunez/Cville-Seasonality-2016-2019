## declare job id
args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])
### ---> launch array from 1-101

### libraries
library(data.table)
library(foreach)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)

### open GDS
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


### samps
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]




### load glm.out
fl <- list.files("/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/GLM_new_mods_data", full.names=T)

### this load one job set at a time
load(fl[jobId])

gsub("/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/GLM_new_mods_data/VA_mod.newmods.","", fl[jobId]) -> perm1
gsub(".Rdata","", perm1) -> perm_num
perm_id = as.numeric(perm_num)

###############
### RNVP ###
###############


#### Calculate RNPval
out.glm = foreach(mod.i=0:11,
                  .combine = "rbind"
                  )%do%{
  
  tmp.rn <- out %>%
    filter(chr %in% c("2L", "2R", "3L", "3R") ) %>%
    filter(mod == mod.i) %>%
    .[complete.cases(p_lrt),] 
  
  tmp.rn %>%
    arrange(p_lrt) %>%
    mutate(rank = seq(from = 1, to = dim(tmp.rn)[1] )) %>%
    mutate(rnp.clean = rank/dim(tmp.rn)[1])
}

###############
### windows ###
###############
# generate a master index for window analysis
### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(out.glm, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
  tmp <- out.glm %>%
        filter(chr == chr.i)
  
  data.table(chr=chr.i,
             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
}

wins[,i:=1:dim(wins)[1]]

dim(wins)

### run windows
### run windows
### run windows

setkey(out.glm, chr, pos)
head(out.glm)

### start the summarization process
win.out <- foreach(win.i=c(1:dim(wins)[1]), 
                   .errorhandling = "remove",
                   .combine = "rbind"
                   )%dopar%{
                     
  message(paste(win.i, dim(wins)[1], sep=" / "))
                     
                     
  win.tmp <- out.glm[J(data.table(chr=wins[win.i]$chr, 
                              pos=wins[win.i]$start:wins[win.i]$end, 
                              key="chr,pos")), nomatch=0]
  
  
  #### Calculate Z score
  win.tmp[,Z:=qnorm(p_lrt, 0, 1)]
  #### Calculate Z rnp score
  win.tmp[,rnpZ:=qnorm(rnp.clean, 0, 1)]
  
  
  seqSetFilter(genofile, 
               variant.id=unique(win.tmp$variant.id),
               sample.id=samps[locality=="VA_ch"][year>= 2016 ]$sampleId)
  
  #obtain AFs 
  af <- seqGetData(genofile, "annotation/format/FREQ")
  f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), 
                      variant.id=seqGetData(genofile, "variant.id"))
  
  #merge AFs with object
  win.tmp <- merge(win.tmp, f.hat, by="variant.id")
  win.tmp[,het:=2*fhat*(1-fhat)]
  
  #thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))[,1]
  
  thrs <- c(0.05)
  
  tmpo <- foreach(pr.i=thrs, 
                  .errorhandling="remove", 
                  .combine="rbind")%do%{
    win.tmp %>% 
    group_by( label, year, start, end, chr ) %>%
    filter(!is.na(rnp.clean)) %>%
    summarize(
      pos_mean = mean(pos),
      pos_min = min(pos),
      pos_max = max(pos),
      win=win.i,
      pr=pr.i,
      perm=perm_id,
      perm_type=ifelse(perm_id==0, "real","permuted"),
      rnp.pr=c(mean(rnp.clean<=pr.i)),
      rnp.binom.p=c(binom.test(sum(rnp.clean<=pr.i), 
                               length(rnp.clean), pr.i)$p.value),
      #rnChi.pr=c(mean(rnChi.clean<=pr.i)),
      #rnChi.binom.p=c(binom.test(sum(rnChi.clean<=pr.i), 
                                 #length(rnChi.clean), pr.i)$p.value),
      wZa=sum(het*Z)/(sqrt(sum(het^2))),
      wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
      rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
      rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
      min.p.lrt=min(p_lrt),
      #min.rank=min(rp.clean),
      min.rnp=min(rnp.clean),
      nSNPs = n(),
      sum.rnp=sum(rnp.clean<=pr.i),
    ) %>% 
    mutate(
    invName=case_when(
      chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
      chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
      chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
      chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
      chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
      chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
      TRUE ~ "noInv"
    ))
  }
  tmpo
}

### save
out_folder <- "/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/GLM_window_analysis"

message(paste(out_folder, "/VA_ch_window_", unique(out.glm$perm) , ".Rdata", sep=""))

save(win.out, 
     file=paste(out_folder, "/VA_ch_window_", unique(out.glm$perm) , ".Rdata", sep=""))

