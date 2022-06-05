### module load intel/18.0 intelmpi/18.0 R/3.6.3; R
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

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
  library(lme4)

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
    glm.out <- foreach(i=snp.dt[job==jobId]$id, .errorhandling="remove")%dopar%{
      ### counter and testers
        # i<-291713; i<-196289
        message(paste(which(i==snp.dt[job==jobId]$id), snp.dt.ag[job==jobId]$N, sep=" / "))

      ### get allele frequencies
        data <- getData(chr=snp.dt[J(i)]$chr,
                        start=snp.dt[J(i)]$pos,
                        end=snp.dt[J(i)]$pos,
                        samplesUse=samps$sampleId)

      ### permutations
        o <- foreach(perm.i=0:10, .combine="rbind")%do%{
          if(perm.i==0) {
            ord.tmp <- samps[,list(aveTemp= aveTemp,
                                   year=year,
                                   sampleId=sampleId), list(locality)]

          } else {
            set.seed(perm.i)
            ord.tmp <- samps[,list(aveTemp=as.numeric(sample(as.character(aveTemp))),
                                   year=as.numeric(sample(as.character(year))),
                                    sampleId=sampleId), list(locality)]
            #ord.tmp[1:5]
          }

          data.use <- merge(data, ord.tmp, by="sampleId")

          ### foreach UseSet that passes filter
            setkey(snp.dt, id)
            #usePops <- melt(snp.dt[J(i), unique(data.use$locality), with=F], measure.vars=unique(data.use$locality))

            ###
            t0 <-  glmer(af_nEff~1 + (1|locality),
                          data=data.use[af_nEff>0 & af_nEff<1],
                          family=binomial(),
                          weights=data.use[af_nEff>0 & af_nEff<1]$nEff)

            t1 <-  glmer(af_nEff~I(aveTemp/10) + (1|locality),
                          data=data.use[af_nEff>0 & af_nEff<1],
                          family=binomial(),
                          weights=data.use[af_nEff>0 & af_nEff<1]$nEff)
            t1.aveTemp.sum <- summary(t1)
            aveTemp.aov <- anova(t1, t0)

            data.table(variant.id=i, chr=data.use[1]$chr, pos=data.use[1]$pos, locality="all", perm=perm.i,
                        mod="aveTemp", year=NA, delta=30,
                         p=coef(t1.aveTemp.sum)[-1,4], b=coef(t1.aveTemp.sum)[-1,1], a=coef(t1.aveTemp.sum)[1,1],
                         chisq=aveTemp.aov$Chisq[2], AIC=AIC(t1), p.lrt=aveTemp.aov$Pr[2],
                         df=aveTemp.aov$Df[2],
                         nSamps=dim(data.use[!is.na(af_nEff) & af_nEff>0 & af_nEff<1])[1],
                         f.hat=mean(data.use[!is.na(af_nEff) & af_nEff>0 & af_nEff<1]$af_nEff))
        }

      ### return
        o
    }
    glm.out <- rbindlist(glm.out)
    glm.out[perm==0][order(p.lrt)]
    glm.out[perm==0][order(p)]

### split and save
  ### first, make sure that directories are there
  # system("rm -R /scratch/aob2x/dest_glm")
  mainDir="/scratch/aob2x/dest_glmer"
  dir.create(file.path(mainDir), showWarnings = FALSE)

  setkey(glm.out, perm)

  foreach(perm.i=c(0:max(glm.out$perm)))%do%{
    # locale.i<-"VA_ch"; perm.i<-0
    subDir=paste("perm", perm.i, sep="")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    tmp <- glm.out[J(data.table(perm=perm.i))]

    save(tmp,
        file=paste(mainDir, "/", subDir, "/all;", perm.i, ";aveTemp_glmer;", jobId, ".Rdata", sep=""))

  }

  ## system("ls -lh /scratch/aob2x/dest_glm/*")
