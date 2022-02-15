# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1]) - 1

## jobId=0


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(SeqArray)
  #library(metap)
  library(lubridate)

### load glm.out
  fl <- list.files("/scratch/aob2x/dest_glm_morePerms_nested_qb/processedGLM/", "glm.out", full.names=T)
  fl <- fl[grepl(paste("_", jobId, ".Rdata", sep=""), fl)]

  load(fl[grepl("VA_ch", fl)])
  message("loading VA data")
  glm.out.va <- glm.out


  fl <- fl[!grepl("VA_ch", fl)]
  message(fl)

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

### open GDS
  message("opening gds")
 genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


###############
### windows ###
###############

### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  setkey(glm.out.va, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%dopar%{
      tmp <- glm.out.va[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

### iterate across other GLM tests
foreach(fl.i=fl)%do%{
  ### load second window
    # fl.i <- fl[6]
    message(fl.i)
    load(fl.i)

  ### merge
    glm.out.va <- glm.out.va[mod=="aveTemp+year_factor"]
    glm.out <- glm.out[mod=="aveTemp+year_factor"]

    setkey(glm.out.va, chr, pos)
    setkey(glm.out, chr, pos)

    message("merge")
    m <- merge(glm.out.va, glm.out)
    m <- m[!is.na(rnp.clean.x) & !is.na(rnp.clean.y)]
    m[,rnp.clean.x:=rank(p.lrt.x) / (length(p.lrt.x) + 1)]
    m[,rnp.clean.y:=rank(p.lrt.y) / (length(p.lrt.y) + 1)]

  ### run windows
    setkey(m, chr, pos)

    win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
      # win.i <- 200
      message(paste(win.i, dim(wins)[1], sep=" / "))
      win.tmp <- m[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
      foreach(thr=c(0.01, 0.05), .combine="rbind", .errorhandling="remove")%do%{
        tmp <- win.tmp[!is.na(rnp.clean.x),
                  list(TT=sum(rnp.clean.x <= thr & rnp.clean.y <= thr),
                      TF= sum(rnp.clean.x <= thr & rnp.clean.y >= thr),
                      FT= sum(rnp.clean.x >= thr & rnp.clean.y <= thr),
                      FF= sum(rnp.clean.x >= thr & rnp.clean.y >= thr),
                      gr=c("obs"),
                      nSNPs=.N, i=win.i,
                      start.bp=min(win.tmp$pos), stop.bp=max(win.tmp$pos), win.i=win.i, thr=thr,
                      invName=paste(unique(invName.x), collapse=";")),
                  list(mod=mod.x, chr, delta=delta.x, locality.x, locality.y, perm=perm.x)]
        tmp[,fet.p:=fisher.test(matrix(c(TT, TF, FT, FF), nrow=2, byrow=T))$p.value]
        tmp[,or:=fisher.test(matrix(c(TT, TF, FT, FF), nrow=2, byrow=T))$estimate]
        return(tmp)
      }
    }

    win.out <- rbindlist(win.out,fill=T)
    win.out <- merge(win.out, wins, by="i")
    save(win.out, file=paste("/scratch/aob2x/dest_glm_morePerms_nested_qb/join_windowAnalysis/joint_window.dt.",
                            m[1]$locality.y, "_", m[1]$perm.x, ".Rdata", sep=""))

    table(win.out$chr.x, win.out$fet.p<.05)
  }
