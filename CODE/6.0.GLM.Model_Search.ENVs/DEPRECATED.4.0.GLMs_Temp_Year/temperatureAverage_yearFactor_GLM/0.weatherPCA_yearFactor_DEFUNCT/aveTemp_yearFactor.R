### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

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

####################################
### sample info and weather data ###
####################################

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

################
### SNP data ###
################

  ### open GDS for common SNPs (PoolSNP)
   genofile.orig <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
   genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds")

  ### common SNP.dt

    #snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
    #                      pos=seqGetData(genofile, "position"),
    #                      nAlleles=seqGetData(genofile, "$num_allele"),
    #                      id=seqGetData(genofile, "variant.id"))

    #snp.dt <- snp.dt[nAlleles==2]
    #seqSetFilter(genofile, snp.dt$id)
    #snp.dt[,mean.af:=seqGetData(genofile, "annotation/info/AF")$data]
    #snp.dt <- snp.dt[mean.af>.1 & mean.af<.9]

    #setkey(snp.dt, chr)
    #snp.dt <- snp.dt[J(c("2L", "2R", "3L", "3R", "X"))]

    load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_Filter_Metadat.Rdata")
    snp.dt <- snp.dt_metadata; rm(snp.dt_metadata)
    snp.dt <- snp.dt[mean.af>.1 & mean.af<.9]

  ### assign job vector
    tmp <- rep(c(1:nJobs), each=ceiling(dim(snp.dt)[1]/nJobs))
    snp.dt[,job:=tmp[1:dim(snp.dt)[1]]]
    snp.dt.ag <- snp.dt[,.N,job]

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

  ### sequential GLM
    glm.out <- foreach(i=snp.dt[job==jobId]$id, .errorhandling="remove")%dopar%{
      ### counter and testers
        #i<-181400; perm.i=0
        message(paste(which(i==snp.dt[job==jobId]$id), snp.dt.ag[job==jobId]$N, sep=" / "))

      ### get allele frequencies
        data <- getData(chr=snp.dt[id==i]$chr,
                        start=snp.dt[id==i]$pos,
                        end=snp.dt[id==i]$pos,
                        samplesUse=samps[useSet==T]$sampleId)

        data <- merge(data, temps.summary, by="sampleId")
        data <- data[!is.na(mu)]

      ### foreach UseSet
        mods.out <- foreach(useSet.i=unique(data[useSet==T]$locality), .combine="rbind")%do%{
          # useSet.i <- "VA_ch"
          t0 <- glm(af_nEff~1, data=data[delta==30][locality==useSet.i][useSet==T], family=binomial(), weights=data[delta==30][locality==useSet.i][useSet==T]$nEff)
          t1.year <- glm(af_nEff~as.factor(year), data=data[delta==30][locality==useSet.i][useSet==T], family=binomial(), weights=data[delta==30][locality==useSet.i][useSet==T]$nEff)
          t1.year.sum <- (summary(t1.year))
          year.aov <- anova(t0, t1.year, test="Chisq")

          mods.aveTemp <- foreach(delta.i=c(30), .combine="rbind")%do%{
            #delta.i<-30
            t1.aveTemp <- glm(af_nEff~mu, data=data[delta==delta.i][locality==useSet.i][useSet==T], family=binomial(), weights=data[delta==delta.i][locality==useSet.i][useSet==T]$nEff)
            t1.aveTemp.sum <- (summary(t1.aveTemp))
            aveTemp.aov <- anova(t0, t1.aveTemp, test="Chisq")

            data.table(variant.id=i, chr=data[1]$chr, pos=data[1]$pos, locality=useSet.i,
                      mod="aveTemp", year=NA, delta=delta.i,
                       p=coef(t1.aveTemp.sum)[-1,4], b=coef(t1.aveTemp.sum)[-1,1], a=coef(t1.aveTemp.sum)[1,1],
                       chisq=aveTemp.aov$Deviance[2], AIC=AIC(t1.aveTemp),
                       df=aveTemp.aov$Df[2], nSamps=sum(!is.na(data[delta==delta.i][locality==useSet.i][useSet==T]$nEff)))
          }

          glm.out <- rbind(
            mods.aveTemp,
            data.table(variant.id=i, chr=data[1]$chr, pos=data[1]$pos, locality=useSet.i,
                      mod="year_factor", year=paste(tstrsplit(row.names(coef(t1.year.sum))[-1], ")")[[2]], collapse=";"), delta=NA,
                       p=paste(coef(t1.year.sum)[-1,4], collapse=";"), b=paste(coef(t1.year.sum)[-1,1], collapse=";"), a=coef(t1.year.sum)[1,1],
                       chisq=year.aov$Deviance[2], AIC=AIC(t1.year),
                       df=year.aov$Df[2], nSamps=sum(!is.na(data[delta==30][locality==useSet.i][useSet==T]$nEff)))
          )
          glm.out
        }
        mods.out
    }
    glm.out <- rbindlist(glm.out)


### save
  save(glm.out, file=paste("/scratch/aob2x/lasso_dest/aveTemp_yearFactor_glm/aveTemp_MultiSample_yearFactor_glm", jobId, ".Rdata", sep=""))
