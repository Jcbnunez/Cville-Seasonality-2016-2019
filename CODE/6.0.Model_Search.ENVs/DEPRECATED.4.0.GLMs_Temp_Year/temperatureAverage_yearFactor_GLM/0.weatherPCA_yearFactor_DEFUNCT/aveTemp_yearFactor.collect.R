# ijob -A berglandlab -c5 -p largemem --mem=50G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(ggplot2)
  library(patchwork)


### load inversion
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")
  setnames(inv.dt, "stop", "end")
  setkey(inv.dt, start, end)

### collect results
  fl <- list.files("/scratch/aob2x/lasso_dest/aveTemp_yearFactor_glm/", "aveTemp_MultiSample_yearFactor_glm", full.name=T)

  glm.out <- foreach(fl.i=fl)%dopar%{
    print(fl.i)
    #fl.i <- fl[1]
    load(fl.i)
    return(glm.out)
  }
  glm.out <- rbindlist(glm.out)

### LRT summaries
  #glm.out.lrt <- glm.out[chr!="X",
  #                      list(chisq=mean(chisq),
  #                             df=mean(df),
  #                             nSamps=mean(nSamps)),
  #                      list(variant.id, chr, pos, mod, delta, useSet)]

  glm.out.lrt <- glm.out; rm(glm.out)
  glm.out.lrt[,p.lrt:= 1 - pchisq(chisq, df)]

  glm.out.lrt[chr=="2L" & pos>= 2225744 & pos<=13154180, inv:="In(2L)t"]
  glm.out.lrt[is.na(inv), inv:="none"]

  lrt.rank <- glm.out.lrt[,list(r=rank(-chisq), l=length(chisq), variant.id, pa=p.adjust(p.lrt, "BH")),
                list(mod, delta, locality)]

  #table(lrt.rank$mod, lrt.rank$pa<.15, lrt.rank$locality)

  setkey(lrt.rank, variant.id, mod, delta, locality)
  setkey(glm.out.lrt, variant.id, mod, delta, locality)
  glm.out.lrt <- merge(glm.out.lrt, lrt.rank)
  glm.out.lrt[,rnp:=r/(l+1)]

  #table(glm.out.lrt$mod, glm.out.lrt$pa<.05)

  load("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_Filter_Metadat.Rdata")

  setkey(glm.out.lrt, chr, pos)
  setkey(snp.dt_metadata, chr, pos)

  glm.out.lrt <- merge(glm.out.lrt, snp.dt_metadata)

  save(glm.out.lrt, file="~/glm.out.lrt.Rdata")
  #ggplot(data=glm.out.lrt, aes(p.lrt)) + geom_histogram() + facet_grid(mod~chr+inv)
  rm(glm.out)

### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  setkey(snp.dt_metadata, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
      tmp <- snp.dt_metadata[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

### define windows
  registerDoMC(5)
  setkey(glm.out.lrt, chr, pos)
  win.out <- foreach(win.i=c(1:dim(wins)[1]))%dopar%{
    # win.i <- 50
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out.lrt[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

    win.tmp[,list(rnp.pr=c(mean(rnp<=0.05), mean(rnp<=0.01)),
                  pr=c(0.05, .01),
                  mean.chisq=c(mean(chisq)),
                  median.chisq=c(median(chisq)),
                  gr=c("obs"),
                  nSNPs=.N, i=win.i),
              list(mod, chr, inv, delta, locality)]

  }

  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i")
  win.out.exp <- data.table(median.chisq=qchisq(.5, c(1,4)),
                            uci.chisq=qchisq(.025, c(1,4)),
                            lci.chisq=qchisq(.975, c(1,4)),
                            mod=c("aveTemp", "year_factor"))

  win.out.exp

  save(win.out, win.out.exp, file="~/win_out_100K_multiPop.Rdata")

###
  library(ggplot2)
  library(data.table)
  library(patchwork)
  load("~/win_out_100K_multiPop.Rdata")

  pr.plot <-
  ggplot(data=win.out[nSNPs>100][chr.x!="X"][gr=="obs"][pr==.05][delta==30 | is.na(delta)][locality=="VA_ch"]) +
  geom_boxplot(aes(x=interaction(chr.x, inv), y=rnp.pr, fill=as.factor(chr.x), group=interaction(chr.x, inv))) +
  facet_wrap(locality~mod, scales="free_y") +
  geom_hline(yintercept=0.05) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
  xlab("") + ylab("Frac. SNPs with \nRank Normalized P <= 0.05)") + ggtitle("Distribution of signal")


  chisq.plot <-
  ggplot(data=win.out[nSNPs>100][chr.x!="X"][gr=="obs"][pr==.05][delta==30 | is.na(delta)][locality=="VA_ch"]) +
  geom_boxplot(aes(x=interaction(chr.x, inv), y=median.chisq, fill=as.factor(chr.x), group=interaction(chr.x, inv))) +
  facet_wrap(locality~mod, scales="free_y") +
  geom_hline(data=win.out.exp, aes(yintercept=median.chisq)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
  xlab("") + ylab("Median Chisq") + ggtitle("Strength of signal")

  strength.abundance <- chisq.plot + pr.plot  + plot_annotation(tag_levels="A")

  ggsave(strength.abundance, file="~/st_ab.VA_ch.png", h=4, w=8)

### MH-plot
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  inv.dt[,y:=c(.3, .3, .3, .35, .4, .3)]
  inv.dt <- inv.dt[,list(delta=c(1,7,14,28,30), pr=.01), list(chr.x, invName, start, stop, y)]
  inv.dt[,a:=1]
  win.out.p <- win.out[,list(p.binom=binom.test(rnp.pr*nSNPs, nSNPs, .05)$p.value), list(i, mod, delta, pr, locality)]

  setkey(win.out, i, mod, delta, pr, locality)
  setkey(win.out.p, i, mod, delta, pr, locality)

  win.out <- merge(win.out, win.out.p)
  win.out[,a:=1]
  win.out[p.binom>1e-5, a:=.5]
  win.out[, a:=-log10(p.binom)]

  mh <- ggplot(data=win.out[mod=="aveTemp"][pr==.05][locality=="VA_ch"],
        aes(x=(start/2 + end/2), y=rnp.pr, color=as.factor(delta))) +
  geom_vline(data=inv.dt, aes(xintercept=start)) + geom_vline(data=inv.dt, aes(xintercept=stop)) +
  geom_segment(data=inv.dt, aes(x=start, xend=stop, y=y/5, yend=y/5), color="black") +
  geom_point(size=1) +
  facet_grid(.~chr.x, scales="free_x") +
  ylim(0, .15)

  ggsave(mh, file="~/mh.VA_ch.png")

### PR plot in other populations
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  setkey(win.out, chr.x)
  win.inv <- foreach(i=1:6)%do%{

    chr.i <- inv.dt[i]$chr.x
    name.i <- inv.dt[i]$invName
    tmp <- win.out[J(chr.i)]

    tmp[,inv:=name.i]

    tmp[,location:="outside"]

    start.i <- inv.dt[i]$start
    stop.i <- inv.dt[i]$stop

    tmp[start>start.i & end<stop.i, location:="inside"]
    tmp[start>(start.i-5e5) & end<(start.i+5e5), location:="brkpoint"]
    tmp[start>(stop.i-5e5) & end<(stop.i+5e5), location:="brkpoint"]

    table(tmp$location)

    tmp
  }
  win.inv <- rbindlist(win.inv)




  pr.plot <-
  ggplot(data=win.inv[chr.x!="X"][gr=="obs"][pr==.01][delta==30][mod=="aveTemp"][nSNPs>25][locality%in%c("VA_ch", "TR_Yes", "DE_Mun")]) +
  geom_boxplot(aes(x=location, y=rnp.pr, fill=as.factor(chr.x), group=interaction(location, locality))) +
  facet_grid(locality~inv) +
  geom_hline(yintercept=0.01) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
  xlab("") + ylab("Frac. SNPs with \nRank Normalized P <= 0.05)") + ggtitle("Distribution of signal")

  ggsave(pr.plot, file="~/pr_world.aveTemp.png")


  dt <- data.table(x=c(1:100), y=c(1:100)+rnorm(100, 0, 75))


  m0 <- lm(y~1, dt)
  m1 <- lm(y~x, dt)

  anova(m0, m1, test="F")



pr.plot <-
ggplot(data=win.out[nSNPs>100][chr.x!="X"][gr=="obs"][pr==.05][delta==30 | is.na(delta)][mod=="aveTemp"]) +
geom_boxplot(aes(x=location, y=rnp.pr, fill=as.factor(chr.x), group=interaction(chr.x, inv))) +
facet_wrap(locality~inv) +
geom_hline(yintercept=0.05) +
theme_bw() +
ylim(0, .2) +
theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
xlab("") + ylab("Frac. SNPs with \nRank Normalized P <= 0.05)") + ggtitle("Distribution of signal")

ggsave(pr.plot, file="~/pr_world.aveTemp.png")


pr.plot <-
ggplot(data=win.out[nSNPs>100][chr.x!="X"][gr=="obs"][pr==.05][delta==30 | is.na(delta)][mod=="year_factor"]) +
geom_boxplot(aes(x=interaction(chr.x, inv), y=rnp.pr, fill=as.factor(chr.x), group=interaction(chr.x, inv))) +
facet_wrap(locality~.) +
geom_hline(yintercept=0.05) +
theme_bw() +
ylim(0, .2) +
theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
xlab("") + ylab("Frac. SNPs with \nRank Normalized P <= 0.05)") + ggtitle("Distribution of signal")

ggsave(pr.plot, file="~/pr_world.year_factor.png")




### old
ggplot(data=glm.out.lrt[p.lrt<=0.05], aes(y=-log10(p.lrt), x=pos, color=chr)) + geom_point() + facet_grid(mod~chr)

dens.plot <-
ggplot() +
geom_density(data=glm.out.lrt[chisq<=25], aes(chisq), fill="red", alpha=.5) +
geom_density(data=glm.out.lrt[chisq<=25], aes(r.chisq), linetype="dashed", fill="green", alpha=.5) +
facet_grid(mod~chr+inv, scales="free_y")

ggsave(dens.plot, file="~/dens_plot.png")

lrt.thr <- foreach(i=0:25, .combine="rbind")%dopar%{
  # i<-4
  message(i)
  glm.out.lrt[,list(nObs=sum(chisq>=i),
                    nSNPs=length(chisq),
                    logpr=pchisq(i, mean(df), log=T, lower.tail=F),
                    nExp=length(chisq) * 10^pchisq(i, mean(df), log=T, lower.tail=F),
                    chisqThr=i),
              list(chr, mod, inv)]
}

lrt.thr[,en:=(nObs)/nExp]

ggplot(data=lrt.thr, aes(x=chisqThr, y=log2(en), group=mod, color=mod)) +
geom_line() +
facet_grid(~chr+inv)

summary(lm(chisq~chr+inv, glm.out.lrt[mod=="year_factor"]))
