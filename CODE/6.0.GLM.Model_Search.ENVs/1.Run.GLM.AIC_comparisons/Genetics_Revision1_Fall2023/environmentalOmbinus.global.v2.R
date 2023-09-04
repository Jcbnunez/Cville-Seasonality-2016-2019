# ijob -A berglandlab_standard -c4 -p dev --mem=20G

### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])
  nJobs=as.numeric(args[2])

  ## jobId=1609; nJobs=5000

### libraries
  library(data.table)
  #library(gdata)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  #library(glmnet)
  library(doMC)
  registerDoMC(8)
  library(fastglm)

  setwd("/scratch/aob2x")

####################################
### sample info and weather data ###
####################################
  ### weather data
    load("Overwintering_18_19/EnvironmentalOmnibus_Global/nasa_power.weather.mod.Rdata")

  ### samps
    samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
    samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
    samps[nchar(month)==1, month:=paste("0", month, sep="")]
    samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
    samps[nchar(day)==1, day:=paste("0", day, sep="")]
    samps <- samps[set!="dgn"]
    samps[,Date:=date(paste(year, month, day, sep="-"))]

    samps[locality=="UA_od", locality:="UA_Ode"]

    samps <- merge(samps, weather.ave, by="sampleId")

### define five population sets: West-Coast NA; East-Coast NA; Western Europe; Eastern Europe; Cville
### also  clean up based on other flags from Kapun et al
  dest <- fread("DEST_freeze1/populationInfo/samps_10Nov2020.csv")
  samps <- merge(samps, dest[status=="Keep"][propSimNorm<.05][,"sampleId"], by="sampleId", all.x=T)

  cluster <- fread("DEST_freeze1/populationInfo/Cluster_Assingment/DEST_Sample_clusters.txt")
  samps <- merge(samps, cluster[,c("sampleId", "Continental_clusters"), with=F], by="sampleId", all.x=T)

  samps[,cluster:=Continental_clusters]
  samps[long< -110, cluster:="2.North_America_W"]
  samps[long> -110 & Continental_clusters=="2.North_America", cluster:="2.North_America_E"]
  samps[locality=="VA_ch", cluster:="5.Cville"]

  samps[locality%in%c("IL_ur", "KA_to", "MI_bh", "ON_su", "WI_cp"), cluster:="2.North_America_Mid"]


  setkey(samps, sampleId)
  table(samps[!duplicated(samps)]$cluster)

  samps <- samps[cluster%in%c("1.Europe_W", "2.North_America_E", "3.Europe_E", "5.Cville")]

  samps.ag <- na.omit(samps[mod==1][,list(.N), list(cluster, locality)])
  samps.ag[,list(numPops=.N), list(N, cluster)][order(N)]
  samps.ag[,list(numPops=sum(N)), list(cluster)]


################
### SNP data ###
################

  ### open GDS for common SNPs (PoolSNP)
   genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

  ### common SNP.dt
    load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
    setkey(snp.dt, id)

    use <- apply(snp.dt[,c("VA_ch"), with=F],
                  1, any)
    snp.dt <- snp.dt[use]
    setkey(snp.dt, id)

    snp.dt <- snp.dt[N==0 & cm_mb>0 & !is.na(cm_mb) & chr!="X"]

  ### assign job vector
    tmp <- rep(c(1:nJobs), each=ceiling(dim(snp.dt)[1]/nJobs))
    snp.dt[,job:=tmp[1:dim(snp.dt)[1]]]
    snp.dt.ag <- snp.dt[,.N,job]

    snp.dt.ag

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

        if(class(dp)[1]!="SeqVarDataList") {
          dp.list <- list()
          dp.list$data <- dp
          dp <- dp.list
        }

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

    anovaFun <- function(m1, m2) {
      #  m1 <- t0; m2<- t1.year
      ll1 <- as.numeric(logLik(m1))
      ll2 <- as.numeric(logLik(m2))

      parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))

      chisq <- -2*(ll1-ll2)

      1-pchisq(chisq, parameter)

    }

  ### sequential GLM
    # jobId <- 81
    glm.method <- 0

    glm.out <- foreach(i=snp.dt[job==jobId]$id, .errorhandling="remove")%dopar%{
      ### counter and testers
        #i<-310918; jobId=450
        # i<-4903; jobId=1
        #i <- 69; i <- snp.dt[job==jobId][69]$id
        message(paste(which(i==snp.dt[job==jobId]$id), snp.dt.ag[job==jobId]$N, sep=" / "))

      ### get allele frequencies
        setkey(snp.dt, id)

        data <- getData(chr=snp.dt[J(i)]$chr,
                        start=snp.dt[J(i)]$pos,
                        end=snp.dt[J(i)]$pos,
                        samplesUse=unique(samps$sampleId))
        setkey(data, sampleId)
        data <- data[!duplicated(data)]
        #data <- merge(data, weather.ave, by="sampleId")
        samps.ag <- samps[,list(cluster=unique(cluster)), sampleId]

        data <- merge(data, samps.ag, by="sampleId")

      ### for each population cluster
        o <- foreach(cluster.i=na.omit(unique(data$cluster)), .errorhandling="remove")%do%{
          # cluster.i <- "5.Cville"

          ### permutations
          o <- foreach(perm.i=0:100, .combine="rbind")%dopar%{
            message(paste(cluster.i, perm.i, sep=" / "))

            ### first the year only model
              if(perm.i==0) {
                ord.tmp <- samps[cluster==cluster.i,list(year=year,
                                       sampleId=sampleId,
                                       orig.sampleId=sampleId),
                                  list(locality=locality, mod)]
              } else {
                ord.tmp <- foreach(mod.i=unique(samps$mod), .combine="rbind")%do%{
                  set.seed(perm.i)
                  samps[mod==mod.i & cluster==cluster.i,list(year=year, mod=mod.i, locality=locality,
                                                             sampleId=sample(sampleId, replace=F), orig.sampleId=sampleId),
                          ]
                }
              }

              setkey(data, sampleId)
              setkey(ord.tmp, sampleId)
              data.use <- merge(data, ord.tmp)
              #data.use <- data.use[locality=="VA_ch"][!is.na(af_nEff)]

              dl <- melt(data.use[!is.na(af_nEff)], id.vars=c("sampleId", "orig.sampleId", "af_nEff", "nEff", "locality", "mod", "year", "cluster"))
              dl[,year_pop:=paste(year, locality, sep=":")]

              tmp <- dl[mod==1]
              y <- tmp$af_nEff
              X.null <- model.matrix(~1, tmp)
              X.year <- model.matrix(~as.factor(year_pop), tmp)
              t0 <- fastglm(x=X.null, y=y, family=binomial(), weights=tmp$nEff, method=glm.method)
              t1.year <- fastglm(x=X.year, y=y, family=binomial(), weights=tmp$nEff, method=glm.method)

              check <- tmp$sampleId

            ### now the environmental models

              data.use <- merge(data,
                                samps[cluster==cluster.i,list(year=year, sampleId=sampleId,orig.sampleId=sampleId,
                                  temp.ave,
                                  temp.var,
                                  temp.max,
                                  temp.min,
                                  temp.propMax,
                                  temp.propmin,
                                  humidity.ave,
                                  humidity.var,
                                  precip.ave,
                                  precip.var),
                                   list(locality=locality, mod)],
                                by="sampleId")


              dl <- melt(data.use[!is.na(af_nEff)], id.vars=c("sampleId", "orig.sampleId", "af_nEff", "nEff", "locality", "mod", "year", "cluster"))
              dl[,year_pop:=paste(year, locality, sep=":")]


              mods.out <- foreach(mod.i=unique(dl$mod), .combine="rbind")%do%{
                foreach(var.i=unique(dl$variable), .combine="rbind")%do%{
                  # mod.i <-4 ; var.i="temp.max"
                  tmp <- dl[mod==mod.i][variable==var.i]
                  setkey(tmp, sampleId)
                  tmp <- tmp[!duplicated(tmp)]
                  tmp[,origValue:=value]
                  if(perm.i!=0) {
                    set.seed(perm.i)
                    tmp[,value:=sample(value)]
                  }

                  X.year.noperm<- model.matrix(~as.factor(year_pop), tmp)
                  t1.year.noperm <- fastglm(x=X.year.noperm, y=y, family=binomial(), weights=tmp$nEff, method=glm.method)

                  X.year.var <- model.matrix(~as.factor(year_pop)+value, tmp)
                  t1.year.var <- fastglm(x=X.year.var, y=y, family=binomial(), weights=tmp$nEff, method=glm.method)

                  data.table(mod=c( mod.i),
                             variable=var.i,
                             AIC=c(AIC(t1.year.var)),
                             b_temp=last(t1.year.var$coef), se_temp=last(t1.year.var$se),
                             nObs=dim(tmp[mod==mod.i])[1],
                             p_lrt=anovaFun(t1.year.noperm, t1.year.var),
                             check=all(tmp$sampleId==check))
                }
              }

              mods.out <- rbind(data.table(mod=c(-1, 0), AIC=c(AIC(t0), AIC(t1.year)),
                                            variable=c("null", "pop_year"),
                                            b_temp=c(NA, NA), se_temp=c(NA,NA),
                                            nObs=length(y),
                                            p_lrt=c(NA, anovaFun(t0, t1.year)),
                                            check=NA),
                                mods.out)

              mods.out[,variant.id:=i]
              mods.out[,perm:=perm.i]
              mods.out[,cluster:=cluster.i]
              return(mods.out)
          }
          return(o)
        }

      ### return
        rbindlist(o)
    }
    glm.out <- rbindlist(glm.out)

    #glm.out[,list(var=variable[which.min(AIC)], mod=mod[which.min(AIC)]), list(variant.id, perm)][perm==0]

### split and save
  message("saving")

  ### first, make sure that directories are there
  # system("rm -R /scratch/aob2x/dest_glm")
  mainDir="/scratch/aob2x/environmental_ombibus_global_permV2"
  dir.create(file.path(mainDir), showWarnings = FALSE)


  save(glm.out, file=paste(mainDir, "/job", jobId, ".Rdata", sep=""))

  message("done")
  ## system("ls -lh /scratch/aob2x/dest_glm/*")
