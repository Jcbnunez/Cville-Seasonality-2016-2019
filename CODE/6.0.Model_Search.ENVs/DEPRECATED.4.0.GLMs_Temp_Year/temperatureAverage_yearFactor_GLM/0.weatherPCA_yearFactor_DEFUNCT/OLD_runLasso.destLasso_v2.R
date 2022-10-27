### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])
  nJobs=as.numeric(args[2])

  ## jobId=44; nJobs=999

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(glmnet)
  library(doMC)
  registerDoMC(4)

####################################
### sample info and weather data ###
####################################

  ### load weather data
#    load(file="/project/berglandlab/alan/weather_impute.Rdata")
    load(file="~/weather_impute.Rdata")

    wm.impute[,term:=paste(variable, dateDelta, sep="_")]

  ### samps
#    samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
    samps <- fread("~/DEST_10Mar2021_POP_metadata.csv")

    samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
    samps[nchar(month)==1, month:=paste("0", month, sep="")]
    samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
    samps[nchar(day)==1, day:=paste("0", day, sep="")]

    setkey(samps, sampleId)
    samps <- samps[J(unique(wm.impute$sampleId))]

  ### make thermal limit model matrix
    wm.hot <- foreach(hot.thr=seq(from=200 ,to=350, by=10), .combine="rbind")%do%{
      wm.impute[,list(nHot=sum(value[variable=="tmax"]>=hot.thr),
                      nCold=sum(value[variable=="tmax"]<=hot.thr),
                      hot.thr=paste("hotThr_", hot.thr, sep="")),
                  list(sampleId)]
    }

    wm.hot[,hot.thr:=factor(hot.thr, levels=paste("hotThr_", seq(from=-100 ,to=350, by=10), sep=""))]
    setkey(wm.hot, sampleId)
    wm.hot <- merge(wm.hot, samps, by="sampleId")

    ww.cville <-dcast(wm.hot[J(samps[locality=="VA_ch"][year>2012]$sampleId)],
                      sampleId~hot.thr, value.var="nHot")

    year.cville <- samps[locality=="VA_ch"][year>2012][,c("year", "sampleId", "yday")]
    year.cville[,Y2015:=as.numeric(year==2015)*10]
    year.cville[,Y2016:=as.numeric(year==2016)*10]
    year.cville[,Y2017:=as.numeric(year==2017)*10]
    year.cville[,Y2018:=as.numeric(year==2018)*10]

    ww.cville <- merge(ww.cville, year.cville[,-"year",with=F], by="sampleId")

    ggplot(data=wm.hot[grepl("VA_ch", sampleId)],
          aes(x=as.numeric(tstrsplit(hot.thr, "_")[[2]]), y=sampleId, fill=nHot)) + geom_tile()

    ggplot(data=wm.hot[grepl("VA_ch", sampleId)],
          aes(x=as.numeric(tstrsplit(hot.thr, "_")[[2]]), y=as.factor(yday), fill=nHot)) + geom_tile()

    ggplot(data=ww.cville, aes(x=yday, y=hotThr_270)) + geom_point()


    ### reduction
    ww.temps <-dcast(wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"],
                    sampleId~term, value.var="value")

      ww.mat <- as.matrix(ww.temps[,-"sampleId", with=F])
      row.names(ww.mat) <- as.matrix(ww.temps[,"sampleId", with=F])[,1]

      ww.temps.pca <- princomp(ww.mat)
      temps.pca <- data.table(sampleId=row.names(ww.temps.pca$scores),
                              pc1=ww.temps.pca$scores[,1])

      temps.pca <- merge(temps.pca, samps, by="sampleId")

      temps.summary <- wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"][,list(mu=mean(value)), list(sampleId)]

      temps.pca <- merge(temps.pca, temps.summary, by="sampleId")


      pca.yday.plot <- ggplot(data=temps.pca, aes(x=yday, y=pc1, color=year)) + geom_point()

      pca.avaTemp.plot <- ggplot(data=temps.pca, aes(x=mu, y=pc1, color=year)) + geom_point()

      temp.plot <- ggplot(data=wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"],
              aes(x=term, y=sampleId, fill=value)) + geom_tile()


      library(patchwork)
        temp.plot + pca.yday.plot +   pca.avaTemp.plot

################
### SNP data ###
################

  ### open GDS for common SNPs (PoolSNP)
#   genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
    genofile <- seqOpen("~/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

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

  ### assign job vector
    tmp <- rep(c(1:nJobs), each=ceiling(dim(snp.dt)[1]/nJobs))
    snp.dt[,job:=tmp[1:dim(snp.dt)[1]]]

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
        afis

    }

  ### test
    lasso.out <- foreach(i=snp.dt[job==jobId]$id[120:220], .errorhandling="remove")%dopar%{
      ### counter and testers
        message(i)
        #i<-181636; perm.i=0

      ### get allele frequencies
        data <- getData(chr=snp.dt[id==i]$chr,
                        start=snp.dt[id==i]$pos,
                        end=snp.dt[id==i]$pos,
                        samplesUse=samps[locality=="VA_ch"][year>2012]$sampleId)

      ### condition on allele global average allele frequency
        if(mean(data$af_nEff, na.rm=T)>=.1 & mean(data$af_nEff, na.rm=T)<=.9) {
          ### first, run LASSO on training set (Cville)

            m.train <- merge(ww.cville, data, by="sampleId")
            m.train <- m.train[!is.na(af_nEff)]
            m.train <- m.train[af_nEff!=0 & af_nEff!=1]
            setnames(m.train, "yday.y", "yday")
            m.train[,yday:=yday/10]

            ord <- 1:dim(m.train)[1]
            #message(perm.i)

            ### fit training set
            o <- cv.glmnet(y=m.train$af_nEff[ord],
                  x=as.matrix(m.train[,names(ww.cville)[-1], with=F]),
                  famliy=binomial(),
                  weights=m.train$nEff[ord],
                  alpha=1, nfolds=50, grouped=F, standardize = FALSE, relax=F)

            fit <- glmnet(y=m.train$af_nEff[ord],
                  x=as.matrix(m.train[,names(ww.cville)[-1], with=F]),
                  famliy=binomial(),
                  weights=m.train$nEff[ord],
                  alpha=1,
                  lambda=o$lambda.1se, grouped=F, standardize = FALSE, relax=T)


            ### format output
              o.mat <- as.matrix(coef(fit, s = "lambda.1se"))
              term <- rownames(o.mat)
              o.mat <- as.data.table(o.mat)
              o.mat[,term:=term]
              o.mat <- o.mat[term!="(Intercept)"][abs(s1)>0]

              pred <- predict(fit, s=o$lambda.1se, relax=T,
                          newx=as.matrix(m.train[,names(ww.cville)[-1], with=F]))

              o.dt <- data.table(variant.id=i,
                                 af=mean(data$af_nEff, na.rm=T),
                                 perm=perm.i,
                                 lambda=o$lambda.1se,
                                 nSamps.train=dim(m.train)[1],
                                 nTerms.train=dim(o.mat)[1],
                                 terms=paste(o.mat$term, collapse=";"),
                                 coefs=paste(o.mat$s1, collapse=";"),
                                 rho.spearman=cor(qlogis(pred[,1]), qlogis(m.train$af_nEff), method="spearman"),
                                 rho.pearson =cor(qlogis(pred[,1]), qlogis(m.train$af_nEff), method="pearson"))

              return(o.dt)
            }
    }
    lasso.out <- rbindlist(lasso.out)

  ### sequential GLM
    glm.out <- foreach(i=snp.dt[job==jobId]$id[120:320], .errorhandling="remove")%dopar%{
      ### counter and testers
        message(i)
        #i<-181400; perm.i=0

      ### get allele frequencies
        data <- getData(chr=snp.dt[id==i]$chr,
                        start=snp.dt[id==i]$pos,
                        end=snp.dt[id==i]$pos,
                        samplesUse=samps[locality=="VA_ch"][year>2012]$sampleId)

        data <- merge(data, temps.summary, by="sampleId")

        t0 <- glm(af_nEff~1, data=data, family=binomial(), weights=data$nEff)
        t1.year <- glm(af_nEff~as.factor(year), data=data, family=binomial(), weights=data$nEff)
        t1.aveTemp <- glm(af_nEff~mu, data=data, family=binomial(), weights=data$nEff)

        t1.year.sum <- (summary(t1.year))
        t1.aveTemp.sum <- (summary(t1.aveTemp))

        year.aov <- anova(t0, t1.year, test="Chisq")
        aveTemp.aov <- anova(t0, t1.aveTemp, test="Chisq")

        glm.out <- rbind(
          data.table(variant.id=i, chr=data[1]$chr, pos=data[1]$pos,
                    hotThr=NA, mod="aveTemp", year=NA,
                     p=coef(t1.aveTemp.sum)[-1,4], b=coef(t1.aveTemp.sum)[-1,1],
                     chisq=aveTemp.aov$Deviance[2],
                     df=aveTemp.aov$Df[2]),

          data.table(variant.id=i, chr=m.train[1]$chr, pos=m.train[1]$pos,
                    hotThr=NA, mod="year_factor", year=tstrsplit(row.names(coef(t1.year.sum))[-1], ")")[[2]],
                     p=coef(t1.year.sum)[-1,4], b=coef(t1.year.sum)[-1,1],
                     chisq=year.aov$Deviance[2],
                     df=year.aov$Df[2])
          )
        glm.out
    }
    glm.out <- rbindlist(glm.out)




    glm.out.ag <- glm.out[,list(chisq=mean(chisq), df=mean(df)), list(chr, pos, mod, variant.id)]

    glm.out.ag[,pChisq:=1-pchisq(chisq, df)]

    ggplot(data=glm.out.ag, aes(pChisq)) + geom_histogram() + facet_wrap(~mod)
### save


  save(lasso.out, file=paste("/scratch/aob2x/lasso_dest/output_destLasso/outPermCV_StandardizeFalse_Relax", jobId, ".Rdata", sep=""))
