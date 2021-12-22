# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])
  nJobs=as.numeric(args[2])

  ## jobId=44; nJobs=1000

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  #library(glmnet)
  library(doMC)
  registerDoMC(1)


  ### load weather data
    load(file="/project/berglandlab/alan/weather_impute.Rdata")

    wm.impute[,term:=paste(variable, dateDelta, sep="_")]

  ### samps
    samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

    samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
    samps[nchar(month)==1, month:=paste("0", month, sep="")]
    samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
    samps[nchar(day)==1, day:=paste("0", day, sep="")]

    setkey(samps, sampleId)
    samps <- samps[J(unique(wm.impute$sampleId))]

    tab <- samps[,list(.N), list(locality, year)]
    tab.nYear <- tab[N>1, list(nYear=length(unique(year))), list(locality)]

    setkey(samps, locality)
    samps[locality%in%tab.nYear[nYear>1]$locality, useSet:=T]
    samps[is.na(useSet), useSet:=F]
    samps[year<2014 & locality=="VA_ch", useSet:=F]

  ### weather averages
    temps.summary <- wm.impute[,
                              list(mu=c(mean(value[delta<=30]),
                                        mean(value[delta<=28]),
                                        mean(value[delta<=14]),
                                        mean(value[delta<=7]),
                                        mean(value[delta<=1])),
                                   delta=c(30,28,14,7,1)),
                              list(sampleId)]

  ### extract allele frequencies function
    getData <- function(chr="2L", start=14617051, end=14617051, samplesUse=samps$sampleId) {
      # chr="2L"; start=14617051; end=14617051; samplesUse=samps$sampleI

      ### filter to target
        snp.tmp <- data.table(chr=chr, pos=start:end)
        setkey(snp.tmp, chr, pos)
        setkey(snp.dt, chr, pos)
        seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id, sample.id=samplesUse, verbose=T)

      ### get frequencies
        message("Allele Freqs")

        ad <- seqGetData(genofile, "annotation/format/AD")
        dp <- seqGetData(genofile, "annotation/format/DP")

        #if(class(dp)[1]!="SeqVarDataList") {
        #  dp.list <- list()
        #  dp.list$data <- dp
        #  dp <- dp.list
        #}

        af <- data.table(ad=expand.grid(ad$data)[,1],
                         dp=expand.grid(dp$data)[,1],
                         sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                         variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

      ### tack them together
        message("merge")
        #afi <- merge(af, snp.dt1.an, by="variant.id")
        afi <- merge(af, snp.dt, by.x="variant.id", by.y="id")

        afi[,af:=ad/dp]

      ### calculate effective read-depth
        afis <- merge(afi, samps, by="sampleId")

        afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
        afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
        afis[,af_nEff:=round(af*nEff)/nEff]

      ### return
        #afis[,-c("n"), with=F]
        afis

    }

  ################
  ### orig SNP data ###
  ################

    ### open GDS for common SNPs (PoolSNP)
     genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
     #genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds")

    ### common SNP.dt

      snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                            pos=seqGetData(genofile, "position"),
                            nAlleles=seqGetData(genofile, "$num_allele"),
                            id=seqGetData(genofile, "variant.id"))

      snp.dt <- snp.dt[nAlleles==2]
      seqSetFilter(genofile, snp.dt$id)
      snp.dt[,mean.af:=seqGetData(genofile, "annotation/info/AF")$data]
      snp.dt <- snp.dt[mean.af>.1 & mean.af<.9]

      setkey(snp.dt, chr)
      snp.dt <- snp.dt[J(c("2L", "2R", "3L", "3R", "X"))]

    ### get original estimates
      orig.est <- getData(chr="2L", start=4362269, end=4362269)
      orig.est <- merge(orig.est, temps.summary[delta==30], by="sampleId")
      summary(glm(af_nEff~mu, data=orig.est[locality=="VA_ch"][year>2012], family=binomial(), weights=orig.est[locality=="VA_ch"][year>2012]$nEff))

    ### close
      seqClose(genofile)

#### new method
    genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds", allow.duplicate=T)
    load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_Filter_Metadat.Rdata")

    seqResetFilter(genofile)
    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                          pos=seqGetData(genofile, "position"),
                          nAlleles=seqGetData(genofile, "$num_allele"),
                          id=seqGetData(genofile, "variant.id"))

    snp.dt <- snp.dt[nAlleles==2]
    seqSetFilter(genofile, snp.dt$id)
    snp.dt[,mean.af:=seqGetData(genofile, "annotation/info/AF")$data]

    setkey(snp.dt, chr)
    snp.dt <- snp.dt[J(c("2L", "2R", "3L", "3R", "X"))]

    ### get new
    rm.est <- getData(chr="2L", start=4362269, end=4362269)
    rm.est <- merge(rm.est, temps.summary[delta==30], by="sampleId")
    summary(glm(af_nEff~mu, data=rm.est[locality=="VA_ch"][year>2012], family=binomial(), weights=rm.est[locality=="VA_ch"][year>2012]$nEff))


### combine
  com <- merge(orig.est[,c("sampleId", "af_nEff"), with=F], rm.est[,c("sampleId", "af_nEff"), with=F])




cat InterruptedRepeats.bed | awk '{
  if($1=="chr2L" && $2<=4362253 && $3>=4362253) print $0
}'


cat WM_SDust.bed | awk '{
  if($1=="chr2L" && $2<=4362253 && $3>=4362253) print $0
}'



      #load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_Filter_Metadat.Rdata")
      #snp.dt <- snp.dt_metadata; rm(snp.dt_metadata)
      #snp.dt <- snp.dt[mean.af>.1 & mean.af<.9]
