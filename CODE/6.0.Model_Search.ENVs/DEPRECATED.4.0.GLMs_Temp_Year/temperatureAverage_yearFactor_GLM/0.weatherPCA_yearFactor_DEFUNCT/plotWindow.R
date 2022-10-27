# ijob -p standard -A berglandlab_standard -c20 --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  library(tidyr)
  library(SeqArray)
  library(lubridate)

### scp aob2x@rivanna.hpc.virginia.edu:~/windows_wza.Rdata ~/.
  load("~/windows_wza.Rdata")

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

### GLM output
  fl <- c("/project/berglandlab/summarized_dest_glm/glm.out.DE_Bro_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.DE_Mun_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.FI_Aka_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.PA_li_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.TR_Yes_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.UA_Ode_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata"
    )

    glmo <- foreach(fli=fl)%dopar%{
      message(fli)
      load(fli)
      glm.out[mod=="aveTemp"][perm==0]
    }

    glmo <- rbindlist(glmo)
    glmo[,b:=as.numeric(as.character(b))]
    glmo[,p:=as.numeric(as.character(p))]
    glmo[,se:=sqrt(b^2 / qchisq(p, df=1, lower=F))]


### open GDS for common SNPs (PoolSNP)
 genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### common SNP.dt
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

### function
  getData <- function(chr="2L", start=14617051, end=14617051) {
    # chr="2L"; start=14617051; end=14617051

    ### filter to target
      snp.tmp <- data.table(chr=chr, pos=start:end)
      setkey(snp.tmp, chr, pos)
      setkey(snp.dt, chr, pos)
      seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)

    ### get annotations
      message("Annotations")
      tmp <- seqGetData(genofile, "annotation/info/ANN")
      len1 <- tmp$length
      len2 <- tmp$data

      snp.dt1 <- data.table(len=rep(len1, times=len1),
                            ann=len2,
                            id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))

    # Extract data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
      snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                            list(variant.id=id)]

      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
      snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

    ### get frequencies
      message("Allele Freqs")

      ad <- seqGetData(genofile, "annotation/format/AD")
      dp <- seqGetData(genofile, "annotation/format/DP")

      af <- data.table(ad=expand.grid(ad$data)[,1],
                       dp=expand.grid(dp$data)[,1],
                       sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                       variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
      message("merge")
      afi <- merge(af, snp.dt1.an, by="variant.id")
      afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")

      afi[,af:=ad/dp]

    ### calculate effective read-depth
      afis <- merge(afi, samps, by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
      afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      afis[,-c("n"), with=F]
  }

### subset
  annotateData <- function(window=1548)
    win.i <- o[perm==0][mod=="aveTemp"][i%in%window]
    afs <- getData(chr=win.i$chr.x, start=min(win.i$start), end=max(win.i$end))

    ### add in GLM output




















  o[perm==0][mod=="aveTemp"][locality=="UA_Ode"][i%in%c(1548)]
  o[perm==0][mod=="aveTemp"][locality=="TR_Yes"][i%in%c(1548)]


  o[perm==0][mod=="aveTemp"][locality=="TR_Yes"][order(wZa.p)][1:5]
  o[perm==0][mod=="aveTemp"][locality=="UA_Ode"][i%in%c(115,116)]
  o[perm==0][mod=="aveTemp"][locality=="VA_ch"][i%in%c(115,116)]
