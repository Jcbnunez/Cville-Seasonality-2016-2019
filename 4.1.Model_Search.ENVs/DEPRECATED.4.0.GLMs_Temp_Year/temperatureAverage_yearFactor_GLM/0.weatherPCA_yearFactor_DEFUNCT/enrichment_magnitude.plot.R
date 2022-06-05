### libraries
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### load
  load("~/glm.out.lrt.Rdata"

### define windows
  win.bp <- 5e5
  step.bp <- 1e5
  setkey(glm.out.lrt, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
      data.table(chr=chr.i,
                  start=seq(from=min(glm.out.lrt[J(chr.i)]$pos), to=max(glm.out.lrt[J(chr.i)]$pos)-win.bp, by=step.bp),
                  end=seq(from=min(glm.out.lrt[J(chr.i)]$pos), to=max(glm.out.lrt[J(chr.i)]$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

### define windows
  setkey(glm.out.lrt, chr, pos)
  win.out <- foreach(win.i=c(1:dim(wins)[1]))%dopar%{
    # win.i <- 100
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out.lrt[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

    win.tmp[,list(rnp.pr=c(mean(rnp<=0.05)),
                  pr=c(0.05),
                  mean.chisq=c(mean(chisq), mean(rchisq(length(chisq), df[1]))),
                  median.chisq=c(median(chisq), median(rchisq(length(chisq), df[1]))),
                  gr=c("obs", "exp"),
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

###

  pr.plot <-
  ggplot(data=win.out[nSNPs>1000][chr.x!="X"]) +
  geom_boxplot(aes(x=interaction(inv, chr.x), y=rnp.pr, fill=chr.x)) +
  facet_grid(~mod) +
  geom_hline(yintercept=0.05) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
  xlab("") + ylab("Frac. SNPs with \nRank Normalized P <= 0.05)") + ggtitle("Distribution of signal")


  chisq.plot <-
  ggplot(data=win.out[nSNPs>1000][chr.x!="X"][gr=="obs"]) +
  geom_boxplot(aes(x=interaction(inv, chr.x), y=median.chisq, fill=chr.x, group=interaction(inv, chr.x))) +
  facet_wrap(.~mod, scales="free_y", nrow=1) +
  geom_hline(data=win.out.exp, aes(yintercept=median.chisq)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.pos="none") +
  xlab("") + ylab("Median Chisq") + ggtitle("Strength of signal")

  strength.abundance <- chisq.plot + pr.plot  + plot_annotation(tag_levels="A")

ggsave(strength.abundance, file="~/st_ab.png", h=4, w=8)
