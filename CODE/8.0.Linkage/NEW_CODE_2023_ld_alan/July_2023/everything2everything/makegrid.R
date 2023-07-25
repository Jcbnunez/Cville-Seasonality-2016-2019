# ijob -A berglandlab_standard -c20 -p standard --mem=30G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### get command line arguments
  args <- commandArgs(trailingOnly=TRUE)
  jobId <- as.numeric(args[1])
  nJobs <- as.numeric(args[2])

  # nJobs <- 1000; jobId <- 446

### libraries
  library(SNPRelate)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### convert to snprelate GDS
  #vcf.fn <- "/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz"
  gds.fn <- "/scratch/aob2x/M_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.snprelate.gds"

  #snpgdsVCF2GDS(vcf.fn, gds.fn, verbose=T)

### open
  genofile <- snpgdsOpen(gds.fn, allow.fork=TRUE)

### make snp table
  snp.dt <- as.data.table(snpgdsSNPList(genofile, sample.id=NULL))
  snp.dt <- snp.dt[chromosome=="2L"][afreq>.05 & afreq<.95]
  setnames(snp.dt, c("chromosome", "position"), c("chr", "pos"))
  setkey(snp.dt, chr, pos)

### bring in environmental model
  load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata")
  glm.out <- glm.out[perm==0]
  setkey(glm.out, chr, pos)

### merge
  snp.dt <- merge(snp.dt, glm.out, all.x=T)
  table(is.na(snp.dt$perm))

### make windows
  win.size <- 50000
  step.size <- 10000

  wins <- data.table(chr="2L",
                     start=seq(from=min(snp.dt$pos), to=max(snp.dt$pos)-win.size, by=step.size),
                     stop=seq(from=min(snp.dt$pos), to=max(snp.dt$pos)-win.size, by=step.size) + win.size)

  wins[,i:=1:dim(wins)[1]]
  # save(wins, file="~/ld_windows_defn.Rdata")

### make pairs of windows
  wins.pair <- as.data.table(expand.grid(wins$i, wins$i))
  wins.pair <- wins.pair[Var1<=Var2]
  wins.pair[,job:=rep(1:nJobs, each=ceiling(dim(  wins.pair)[1]/nJobs))[1:dim(  wins.pair)[1]]]

### get jobs for this job
  # jobId=50
  wins.pair <- wins.pair[job==jobId]
  setkey(wins.pair, Var1, Var2)
  wins.pair[,x:=1:dim(wins.pair)[1]]
  setkey(wins.pair, x)


  # wins.pair <- wins.pair[1:5]

### iterate
  o <- foreach(i=1:dim(wins.pair)[1])%dopar%{
      # i <- 2754
      win1 <- wins.pair[J(i)]$Var1
      win2 <- wins.pair[J(i)]$Var2


      snp.list1 <- snp.dt[J(data.table(chr=wins[i==win1]$chr, pos=wins[i==win1]$start:wins[i==win1]$stop, key="chr,pos")), nomatch=0]
      snp.list2 <- snp.dt[J(data.table(chr=wins[i==win2]$chr, pos=wins[i==win2]$start:wins[i==win2]$stop, key="chr,pos")), nomatch=0]

      ldmat <- snpgdsLDMat(genofile,
                  sample.id=NULL,
                  snp.id=unique(c(snp.list1$snp.id, snp.list2$snp.id)),
                  slide=0,
                  method=c("corr"),
                  num.thread=1, with.id=TRUE, verbose=TRUE)



      ld.dt <- data.table(r=expand.grid(ldmat$LD)$Var1,
                          use=expand.grid(lower.tri(ldmat$LD))$Var1,
                          snp.id1=rep(ldmat$snp.id, each=length(ldmat$snp.id)),
                          snp.id2=rep(ldmat$snp.id, times=length(ldmat$snp.id)))

      ld.dt <- ld.dt[use==T][snp.id1%in%snp.list1$snp.id][snp.id2%in%snp.list2$snp.id]

      ld.dt <- merge(ld.dt, snp.dt[,c("snp.id", "rnp"),with=F], by.x="snp.id1", by.y="snp.id")
      ld.dt <- merge(ld.dt, snp.dt[,c("snp.id", "rnp"),with=F], by.x="snp.id2", by.y="snp.id")
      ld.dt[,poolOnly:=!is.na(rnp.x) & !is.na(rnp.y)]
      table(ld.dt$genoOnly)

      ld.ag <- foreach(rnp.i=c(.05, .5, 1), .combine="rbind")%do%{
        ld.dt[rnp.x<=rnp.i & rnp.y<=rnp.i, list(meanR2=mean(r^2), medR2=median(r^2), maxR2=max(r^2), prop60=mean(r^2 > .6),
                                                rnp.thr=rnp.i, nMarkers=length(unique(c(snp.id1,snp.id2))), .N), list(poolOnly)]
      }
      ld.ag <- rbind(ld.dt[, list(meanR2=mean(r^2), medR2=median(r^2), maxR2=max(r^2), prop60=mean(r^2 > .6),
                                  rnp.thr=rnp.i, nMarkers=length(unique(c(snp.id1,snp.id2))), .N, poolOnly=F), ],
                     ld.ag, fill=T)

      ld.ag[,nM:=.5+sqrt(1+8*N)/2]
      ld.ag[,win1:=win1]
      ld.ag[,win2:=win2]
      ld.ag[,job:=jobId]
      ld.ag[,start1:=wins[win1]$start]
      ld.ag[,stop1:= wins[win1]$stop]
      ld.ag[,mid1:=  wins[win1]$start/2 +
                     wins[win1]$stop/2]

      ld.ag[,start2:=wins[win2]$start]
      ld.ag[,stop2:= wins[win2]$stop]
      ld.ag[,mid2:=  wins[win2]$start/2 +
                     wins[win2]$stop/2]
      ld.ag[,win.size:=win.size]
      ld.ag[,step.size:=step.size]
      ld.ag
    }

  o <- rbindlist(o)

### save
  save(o, file=paste("/scratch/aob2x/everythingLD/bits/job", jobId, "_", win.size, "_", step.size, ".Rdata", sep=""))
