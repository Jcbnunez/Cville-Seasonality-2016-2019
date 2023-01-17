### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(tidyr)
  library(doMC)
  registerDoMC(5)
#  library(patchwork)
#  library(ggplot2)

### open GDS
  genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### load BED files
  beds.fn <- list.files("~/Overwintering_18_19/RepeatFilterFiles/", ".bed", full.names=T)
  beds.fn <- beds.fn[!grepl("combine", beds.fn)]

  rep.bed <- foreach(i=beds.fn, .combine="rbind")%do%{
    # i<- beds.fn[5]
    message(i)
    tmp <- fread(i, skip=1)
    tmp[,chr:=gsub("chr", "", V1)]
    setnames(tmp, c("V2", "V3"), c("start", "end"))
      ### tmp[chr=="2R"][start<=6481518][end>=6481518]
    tmp <- tmp[,c("chr", "start", "end"), with=F]
    tmp[,repLib:=tstrsplit(i, "/") %>% last %>% gsub(".bed", "", .)]
    tmp # tmp[chr=="2R"][start<=6481518][end>=6481518]
  }
  rep.bed[chr=="2R"][start<=6481518][end>=6481518]

### load basic SNP table
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        id=seqGetData(genofile, "variant.id"),
                        DP=seqGetData(genofile, "annotation/info/DP"))


  snp.dt <- snp.dt[nAlleles==2]
  setkey(snp.dt)
  snp.dt[,start:=pos]
  snp.dt[,end:=pos]

### merge repetitive regions from BED files with SNP table
  setkey(snp.dt, start, end)
  setkey(rep.bed, start, end)

  snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
    #chr.i<-"2L"
    print(chr.i)
    tmp <- foverlaps(snp.dt[chr==chr.i], rep.bed[chr==chr.i])
    tmp.ag <- tmp[,list(.N, libs=paste(repLib, collapse=";"), DP=mean(DP)), list(chr=i.chr, pos, id)]
    tmp.ag[libs=="NA", libs:=NA]
    tmp.ag[is.na(libs), N:=0]
    tmp.ag
  }
  snp.dt <- rbindlist(snp.dt)
  snp.dt[,start:=pos]
  snp.dt[,end:=pos]
  setkey(snp.dt, start, end)

### load in recombination rate file
  recRate <- fread("~/Overwintering_18_19/RecombinationRate/RecRates-All-Chromosomes-100kb.r6.bed")
  setnames(recRate, names(recRate), c("chr", "start", "end", "cm_mb"))
  recRate[,chr:=gsub("chr", "", chr)]
  recRate
  setkey(recRate, start, end)

  snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
    #chr.i<-"2L"
    print(chr.i)
    tmp <- foverlaps(snp.dt[chr==chr.i], recRate[chr==chr.i])
    tmp[,c("chr", "pos", "id", "N", "libs", "cm_mb"), with=F]
  }
  snp.dt <- rbindlist(snp.dt)
  snp.dt[,start:=pos]
  snp.dt[,end:=pos]
  setkey(snp.dt, start, end)

### add in Inversion identifier
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")
  setnames(inv.dt, "stop", "end")
  setkey(inv.dt, start, end)

  snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
    #chr.i<-"3R"
    print(chr.i)
    tmp <- foverlaps(snp.dt[chr==chr.i], inv.dt[chr==chr.i])
    tmp <- tmp[,c("i.chr", "pos", "id", "N", "libs", "cm_mb", "invName"), with=F]
    dim(tmp)
    tmp2 <- tmp[,list(invName=paste(invName, collapse=";")), list(chr=i.chr, pos, id, N, libs, cm_mb)]
    tmp2[invName=="NA", invName:="none"]

    tmp2
  }
  snp.dt <- rbindlist(snp.dt)

### create locality specific filter sets: fix == 0% (no fixed sites), mean AF > 5%, and missing rate <=25%
  samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                        "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                        "VA_ch"=list(locality="VA_ch",     minYear=2014),
                        "PA_li"=list(locality="PA_li",     minYear=2000),
                        "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                        "UA_Ode"=list(locality=c("UA_Ode", "UA_od"), minYear=2010),
                        "TR_Yes"=list(locality="TR_Yes",   minYear=2000))


  setkey(samps, locality)
  popFilter <- foreach(i=1:length(targetLocales))%do%{
    #i <- 1
    #samps[J(targetLocales[[i]]$locality)][year>=targetLocales[[i]]$minYear]

    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt$id, sample.id=samps[J(targetLocales[[i]]$locality)][year>=targetLocales[[i]]$minYear]$sampleId)

    afs <- seqGetData(genofile, "annotation/format/FREQ")

    tmp.dt <- data.table(id=seqGetData(genofile, "variant.id"),
                         meanAF=colMeans(afs$data, na.rm=T),
                         missing=seqMissing(genofile),
                         propFixed=colMeans(afs$data==0, na.rm=T) + colMeans(afs$data==1, na.rm=T))

    tmp.dt[,pass:=F]
    tmp.dt[(meanAF>=0.05 & meanAF<=.95) & propFixed==0 & missing<=0.25, pass:=T]

    #table(tmp.dt$pass)
    tmp.dt[,sampleId:=names(targetLocales)[i]]
    setkey(tmp.dt, id)
    tmp.dt[,c("id", "pass", "sampleId"), with=F]
  }
  popFilter <- rbindlist(popFilter)

  popFilter <- dcast(popFilter, id~sampleId, value.var="pass")

  snp.dt <- merge(snp.dt, popFilter, by="id")

  save(snp.dt, file="~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
