# ijob -A berglandlab_standard -c20 -p standard --mem=30G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### load
  setwd("/scratch/aob2x/ld")
  ld.dat <- fread("plink.ld")
  setkey(ld.dat, CHR_B)

  ld.dat <- ld.dat[J("2L")]
  ld.dat[BP_A<10e6, break_point:="left"]
  ld.dat[BP_A>10e6, break_point:="right"]

### bring in environmental model
  load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata")
  glm.out <- glm.out[perm==0]

### merge
  ld.dat[,pos:=BP_B]
  ld.dat[,chr:=CHR_B]
  setkey(ld.dat, chr, pos)
  setkey(glm.out, chr, pos)

  ld.dat <- merge(ld.dat, glm.out)

  mat <- table(ld.dat$rnp<.0005, ld.dat$R2>.5)
  fisher.test(mat)


### define windows
  win.size <- 25000
  step.size <- 1000

  wins <- data.table(start=seq(from=min(ld.dat$BP_B), to=max(ld.dat$BP_B)-win.size, by=step.size),
                      stop=seq(from=min(ld.dat$BP_B), to=max(ld.dat$BP_B)-win.size, by=step.size) + win.size)

### loop
  setkey(ld.dat, BP_B)
  ld.ag <- foreach(i=1:dim(wins)[1], .combine="rbind")%dopar%{
    message(paste(i, dim(wins)[1], sep=" / "))
    tmp <- ld.dat[J(wins[i]$start:wins[i]$stop)]
    tmp <- tmp[!is.na(BP_A)]
    out.all <- tmp[,list(maf=0, rnp=1, i=i, start=wins[i]$start, stop=wins[i]$stop, r2_mean=mean(R2, na.rm=T), r2_median=median(R2, na.rm=T), r2_max=max(R2, na.rm=T), .N),
        list(break_point)]

    out.nolow <- tmp[MAF_B>.05 & MAF_B<.95,
        list(maf=.05, rnp=1, i=i, start=wins[i]$start, stop=wins[i]$stop, r2_mean=mean(R2, na.rm=T), r2_median=median(R2, na.rm=T), r2_max=max(R2, na.rm=T), .N),
        list(break_point)]


    out.rnp <- foreach(j=c(.0001, .001, .01, .05), .combine="rbind", .errorhandling="remove")%do%{
      out.rnp <- tmp[MAF_B>.05 & MAF_B<.95 & rnp<=j,
          list(maf=.05, rnp=j, i=i, start=wins[i]$start, stop=wins[i]$stop, r2_mean=mean(R2, na.rm=T), r2_median=median(R2, na.rm=T), r2_max=max(R2, na.rm=T), .N),
          list(break_point)]

    }

    rbind(out.all, out.nolow, out.rnp)
  }

### save
  save(ld.ag, file="~/ld_ag.Rdata")


### download and plot
  system("scp aob2x@rivanna.hpc.virginia.edu:~/ld_ag.Rdata ~/.")

library(ggplot2)
library(data.table)

load("~/ld_ag.Rdata")

ld.ag[start<=13193643 & stop>=13193643][order(r2_mean)]

ggplot(data=ld.ag[N>0], aes(x=I(start/2 + stop/2), y=r2_mean), alpha=.05) + geom_line() +
facet_grid(maf + rnp ~break_point, scales="free_y")
