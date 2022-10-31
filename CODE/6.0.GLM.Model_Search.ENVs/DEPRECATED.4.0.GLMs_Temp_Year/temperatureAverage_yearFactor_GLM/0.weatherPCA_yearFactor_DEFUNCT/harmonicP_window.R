# ijob -p standard -A berglandlab_standard -c20 --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(tidyr)
  library(lubridate)
  library(harmonicmeanp)

### load glm output
  load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

### make windows
  glm.out <- glm.out[!is.na(rnp.clean)]
  setkey(glm.out, chr, pos)
  glm.out[,index:=c(1:dim(glm.out)[1])]
  win.bp <- 500
  step.bp <- 250
  setkey(glm.out, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
      tmp <- glm.out[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$index), to=max(tmp$index)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$index), to=max(tmp$index)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

### run windows
  L <- sum(!is.na(glm.out[mod=="aveTemp"]$rnp.clean))
  setkey(glm.out, index)
  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    # win.i <- 126
    message(paste(win.i, dim(wins)[1], sep=" / "))
    #win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
    win.tmp <- glm.out[J(data.table(index=wins[win.i]$start:wins[win.i]$end, key="index")), nomatch=0]


    win.tmp[!is.na(rnp.clean),
              list(rnp.pr=c(mean(rnp.clean<=0.05)),
                  pr=c(0.05),
                  harmonicP=p.hmp(p=na.omit(p.lrt), L=L),
                  stat=hmp.stat(p=na.omit(p.lrt)),

                  mean.chisq=c(mean(chisq)),
                  median.chisq=c(median(chisq)),
                  gr=c("obs"),
                  nSNPs=.N, i=win.i,
                  start.bp=min(win.tmp$pos), stop.bp=max(win.tmp$pos)),
              list(mod, chr, invName, delta, locality)]

  }
  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i")
  win.out[mod=="aveTemp"][order(stat)]


save(win.out, file="~/hp.Rdata")




load("~/hp.Rdata")
library(ggplot2)
library(data.table)


ggplot(win.out[mod=="aveTemp"][nSNPs>200], aes(x=start, y=-log10(stat), color=chr.x)) + geom_line() + facet_grid(~chr.x)
