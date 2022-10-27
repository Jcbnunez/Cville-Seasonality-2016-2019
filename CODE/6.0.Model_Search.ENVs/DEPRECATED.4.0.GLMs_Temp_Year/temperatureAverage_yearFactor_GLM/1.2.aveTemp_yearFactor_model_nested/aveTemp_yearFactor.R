### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])
  nJobs=as.numeric(args[2])

  ## jobId=1; nJobs=1000

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  #library(glmnet)
  library(doMC)
  registerDoMC(5)

####################################
### sample info and weather data ###
####################################
  ### weather data
    load("~/Overwintering_18_19/temperatureAverage_yearFactor_GLM/weatherAve.Rdata")
    setnames(weather.ave, "V1", "sampleId")

  ### samps
    samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
    samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
    samps[nchar(month)==1, month:=paste("0", month, sep="")]
    samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
    samps[nchar(day)==1, day:=paste("0", day, sep="")]
    samps <- samps[set!="dgn"]
    samps[,Date:=date(paste(year, month, day, sep="-"))]

    samps[locality=="UA_od", locality:="UA_Ode"]

    samps <- merge(samps, weather.ave[,c("aveTemp", "sampleId")], by="sampleId")


    targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                          "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                          "VA_ch"=list(locality="VA_ch",     minYear=2014),
                          "PA_li"=list(locality="PA_li",     minYear=2000),
                          "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                          "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                          "TR_Yes"=list(locality="TR_Yes",   minYear=2000))

################
### SNP data ###
################

  ### open GDS for common SNPs (PoolSNP)
   genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

  ### common SNP.dt
    load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
    setkey(snp.dt, id)

    use <- apply(snp.dt[,c("DE_Bro", "DE_Mun", "FI_Aka", "PA_li", "TR_Yes", "UA_Ode", "VA_ch"), with=F],
                  1, any)
    snp.dt <- snp.dt[use]
    setkey(snp.dt, id)

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
        afis[,c("sampleId", "af_nEff", "nEff"), with=F]

    }

  ### sequential GLM
    glm.out <- foreach(i=snp.dt[job==jobId]$id)%dopar%{
      ### counter and testers
        #i<-4
        message(paste(which(i==snp.dt[job==jobId]$id), snp.dt.ag[job==jobId]$N, sep=" / "))

      ### get allele frequencies
        setkey(snp.dt, id)

        data <- getData(chr=snp.dt[J(i)]$chr,
                        start=snp.dt[J(i)]$pos,
                        end=snp.dt[J(i)]$pos,
                        samplesUse=samps$sampleId)

      ### permutations
        o <- foreach(perm.i=0:100, .combine="rbind")%dopar%{
          if(perm.i==0) {
            ord.tmp <- samps[,list(aveTemp= aveTemp,
                                   year=year,
                                   sampleId=sampleId), list(locality)]

          } else {
            set.seed(perm.i)
            ord.tmp <- samps[,list(aveTemp=aveTemp,
                                   year=year,
                                   sampleId=sample(sampleId)),
                              list(locality)]
            #ord.tmp[1:5]
          }

          data.use <- merge(data, ord.tmp, by="sampleId")

          ### foreach UseSet that passes filter
            setkey(snp.dt, id)
            usePops <- melt(snp.dt[J(i), unique(data.use$locality), with=F], measure.vars=unique(data.use$locality))

            mods.out <- foreach(useSet.i=usePops[value==T]$variable, .combine="rbind")%do%{
              # useSet.i <- "VA_ch"
              t0 <- glm(af_nEff~1, data=data.use[locality==useSet.i], family=binomial(), weights=data.use[locality==useSet.i]$nEff)

              t1.year <- glm(af_nEff~as.factor(year), data=data.use[locality==useSet.i], family=quasibinomial(), weights=data.use[locality==useSet.i]$nEff)
              t1.year.sum <- (summary(t1.year))
              year.aov <- anova(t0, t1.year, test="Chisq")

              # t1.aveTemp <- glm(af_nEff~(aveTemp), data=data.use[locality==useSet.i], family=binomial(), weights=data.use[locality==useSet.i]$nEff)
              # t1.aveTemp.sum <- (summary(t1.aveTemp))
              # aveTemp.aov <- anova(t0, t1.aveTemp, test="Chisq")

               t1.aveTemp <- glm(af_nEff~ as.factor(year) + (aveTemp), data=data.use[locality==useSet.i], family=binomial(), weights=data.use[locality==useSet.i]$nEff)
               t1.aveTemp.sum <- (summary(t1.aveTemp))
               aveTemp.aov <- anova(t1.year, t1.aveTemp, test="Chisq")


              rbind(
                data.table(variant.id=i, chr=data.use[1]$chr, pos=data.use[1]$pos, locality=useSet.i, perm=perm.i,
                          mod="year_factor", year=paste(tstrsplit(row.names(coef(t1.year.sum))[-1], ")")[[2]], collapse=";"), delta=30,
                           p=paste(coef(t1.year.sum)[-1,4], collapse=";"), b=paste(coef(t1.year.sum)[-1,1], collapse=";"), a=coef(t1.year.sum)[1,1],
                           chisq=year.aov$Deviance[2], AIC=AIC(t1.year),
                           df=year.aov$Df[2], nSamps=sum(!is.na(data.use[locality==useSet.i]$nEff))),

                data.table(variant.id=i, chr=data.use[1]$chr, pos=data.use[1]$pos, locality=useSet.i, perm=perm.i,
                          mod="aveTemp+year_factor", year=NA, delta=30,
                           p=paste(coef(t1.aveTemp.sum)[-1,4], collapse=";"), b=paste(coef(t1.aveTemp.sum)[-1,1], collapse=";"), a=coef(t1.aveTemp.sum)[1,1],
                           chisq=aveTemp.aov$Deviance[2], AIC=AIC(t1.aveTemp),
                           df=aveTemp.aov$Df[2], nSamps=sum(!is.na(data.use[locality==useSet.i]$nEff)))
              )
            }
            return(mods.out)
        }

      ### return
        o
    }
    glm.out <- rbindlist(glm.out)
    dim(glm.out)

    head(glm.out)

### split and save
  message("saving")

  ### first, make sure that directories are there
  # system("rm -R /scratch/aob2x/dest_glm")
  mainDir="/scratch/aob2x/dest_glm_morePerms_nested"
  dir.create(file.path(mainDir), showWarnings = FALSE)

  setkey(glm.out, locality)

  foreach(locale.i=unique(samps$locality))%do%{
    # locale.i<-"VA_ch"; perm.i<-0
    subDir=locale.i
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    tmp <- glm.out[J(data.table(locality=locale.i))]

    write.csv(tmp,
        file=paste(mainDir, "/", subDir, "/", locale.i, ";allPerm;aveTemp_yearFactor_nested;", jobId, ".csv", sep=""), quote=F, row.names=F)
  }

  message("done")
  ## system("ls -lh /scratch/aob2x/dest_glm/*")
