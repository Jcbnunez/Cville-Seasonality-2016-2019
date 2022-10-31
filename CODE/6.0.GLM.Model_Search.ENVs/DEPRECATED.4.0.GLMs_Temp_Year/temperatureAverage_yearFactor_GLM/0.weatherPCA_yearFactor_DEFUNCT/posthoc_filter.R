

### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)
  library(foreach)


### merge with glm.out
  load("~/Overwintering_18_19/SNP_filtering/snp_dt.Rdata")
  load("~/glm.out.lrt.Rdata")
  setkey(glm.out.lrt, chr, pos)
  setkey(snp.dt, chr, pos)
  glm.out.lrt <- merge(glm.out.lrt, snp.dt)

### rerank after removing bad sites
  lrt.rank <- glm.out.lrt[N==0 & cm_mb>0 & !is.na(cm_mb) & VA_ch==T, list(r.clean=rank(-chisq), l.clean=length(chisq), variant.id),
                list(mod)]

  setkey(lrt.rank, variant.id, mod)
  setkey(glm.out.lrt, variant.id, mod)

  glm.out.lrt <- merge(glm.out.lrt, lrt.rank, all.x=T)
  glm.out.lrt[,rnp.clean:=r.clean/(l.clean+1)]

  #save(glm.out.lrt, file="~/glm_out_lrt.rep.rec.Rdata")

####################################
### some basic tests and figures ###
####################################

  ### temperature test
    tab <- table(sig=glm.out.lrt[mod=="aveTemp"][VA_ch==T]$rnp<0.01,
                rep=glm.out.lrt[mod=="aveTemp"][VA_ch==T]$N>0)
    temp.rep <- fisher.test(tab)
    tab <- table(sig=glm.out.lrt[mod=="aveTemp"][VA_ch==T]$rnp<0.01,
                rec=glm.out.lrt[mod=="aveTemp"][VA_ch==T]$cm_mb==0)
    temp.rec <- fisher.test(tab)

  ### year test
    tab <- table(sig=glm.out.lrt[mod=="year_factor"][VA_ch==T]$rnp<0.01, rep=glm.out.lrt[mod=="year_factor"][VA_ch==T]$N>0)
    year.rep <- fisher.test(tab)
    tab <- table(sig=glm.out.lrt[mod=="year_factor"][VA_ch==T]$rnp<0.01, rep=glm.out.lrt[mod=="year_factor"][VA_ch==T]$cm_mb==0)
    year.rec <- fisher.test(tab)

  ### all inversion, all tests
    fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind")%dopar%{
      # inv.i="2Lt"
      foreach(mod.i=c("aveTemp", "year_factor"), .combine="rbind")%do%{
        print(paste(inv.i, mod.i, sep=" / "))
        # mod.i="aveTemp"
        tab <- table(sig=glm.out.lrt[mod==mod.i][N==0][cm_mb>0 & !is.na(cm_mb)][VA_ch==T]$rnp.clean<0.01,
                    Inv=grepl(inv.i, glm.out.lrt[mod==mod.i][N==0][cm_mb>0 & !is.na(cm_mb)][VA_ch==T]$invName))


        temp.inv <- fisher.test(tab)
        data.table(model=mod.i, factor=inv.i,
                  or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2])
      }
    }

  ### merge
    fet.dt <- rbind(
              fet.invs,
              data.table(model=c(rep("aveTemp", 2), rep("year_factor", 2)),
                         factor=rep(c("Repetitive", "LowRec"), 2),
                          or=c(temp.rep$estimate, temp.rec$estimate,
                               year.rep$estimate, year.rec$estimate),
                          or.lci=c(temp.rep$conf.int[1], temp.rec$conf.int[1],
                                  year.rep$conf.int[1], year.rec$conf.int[1]),
                          or.uci=c(temp.rep$conf.int[2], temp.rec$conf.int[2],
                                   year.rep$conf.int[2], year.rec$conf.int[2]))
              )

    save(fet.dt, file="~/fet.dt")

  ### plot
    load("~/fet.dt")
    ggplot(data=fet.dt, aes(x=factor, y=(or), color=model, group=model)) +
    geom_point(size=4, stat="identity", position=position_dodge(width = .5)) +
    geom_errorbar(aes(ymin=(or.lci), ymax=(or.uci)), width=.1, position=position_dodge(width = .5)) +
    theme_bw()


  ### basic MH plot plot
    inv.dt[,y:=c(15, 15, 15, 15.5, 16, 15)]

    mh.raw <- ggplot(data=glm.out.lrt[p.lrt<.01],
          aes(x=pos, y=-log10(p.lrt))) +
          geom_point(size=1) +
          facet_grid(mod~chr, scales="free_x") +
          ylim(0,17) +
          geom_segment(data=inv.dt, aes(x=start, xend=stop, y=y, yend=y), color="black") +
          ggtitle("No filter")

    mh.filter <- ggplot(data=glm.out.lrt[p.lrt<.01][N==0 & cm_mb>0 & !is.na(cm_mb)],
          aes(x=pos, y=-log10(p.lrt))) +
          geom_point(size=1) +
          facet_grid(mod~chr, scales="free_x") +
          ylim(0,17) +
          geom_segment(data=inv.dt, aes(x=start, xend=stop, y=y, yend=y), color="black") +
          ggtitle("Remove Rep, Remove 0 cm/Mb, remove NA cm/Mb")

    mega <- mh.raw + mh.filter

    ggsave(mega, file="~/mh.full.png", h=6, w=12)


### window plot
  win.bp <- 5e5
  step.bp <- 1e5
  setkey(glm.out.lrt, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
      tmp <- glm.out.lrt[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)


  registerDoMC(5)
  setkey(glm.out.lrt, chr, pos)

  win.out <- foreach(win.i=c(1:dim(wins)[1]))%dopar%{
    # win.i <- 50
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out.lrt[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

    win.tmp[!is.na(rnp.clean),list(rnp.pr=c(mean(rnp.clean<=0.05), mean(rnp.clean<=0.01)),
                  pr=c(0.05, .01),
                  mean.chisq=c(mean(chisq)),
                  median.chisq=c(median(chisq)),
                  gr=c("obs"),
                  nSNPs=.N, i=win.i),
              list(mod, chr, inv)]

  }

  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i")
  win.out.exp <- data.table(median.chisq=qchisq(.5, c(1,4)),
                            uci.chisq=qchisq(.025, c(1,4)),
                            lci.chisq=qchisq(.975, c(1,4)),
                            mod=c("aveTemp", "year_factor"))

  win.out.exp

  setnames(win.out, "chr.x", "chr")

  mh <- ggplot(data=win.out[pr==.01][nSNPs>50],
        aes(x=(start/2 + end/2), y=rnp.pr)) +
  geom_hline(aes(yintercept=.01)) +
  geom_segment(data=inv.dt, aes(x=start, xend=end, y=y/100, yend=y/100), color="black") +
  geom_line(size=1) +
  facet_grid(mod~chr, scales="free_x") +
  ylim(0, .25)

  ggsave(mh, file="~/mh.png")










  save(win.out, win.out.exp, file="~/win_out_500K.filter.Rdata")










mh <- ggplot(data=glm.out.lrt[p.lrt<.25],
      aes(x=pos, y=-log10(p.lrt))) +
#geom_vline(data=inv.dt, aes(xintercept=start)) + geom_vline(data=inv.dt, aes(xintercept=stop)) +
#geom_segment(data=inv.dt, aes(x=start, xend=stop, y=y/5, yend=y/5), color="black") +
geom_point(size=1) +
facet_grid(mod~chr, scales="free_x")

ggsave(mh, file="~/mh.png")








  tab <- table(sig=m$rnp<.05, rep=grepl("SimpleRepeats", m$libs))
  fisher.test(tab)

  tab <- table(sig=m$rnp<.05, rep=grepl("RepeatMasker", m$libs))
  fisher.test(tab)

  tab <- table(sig=m$rnp<.05, rep=grepl("WM_SDust", m$libs))
  fisher.test(tab)

  tab <- table(sig=m$rnp<.05, rep=grepl("InterruptedRepeats", m$libs))
  fisher.test(tab)

  tab <- table(sig=m$rnp<.05, rep=grepl("InterruptedRepeats", m$libs))
  fisher.test(tab)



### year
glmy <- glm.out.lrt[mod=="year_factor"]

setkey(glmy, chr, pos)
my<- merge(snp.dt, glmy)

tab <- table(sig=my$r<100, rep=my$N>0)
fisher.test(tab)

my[r<10]
