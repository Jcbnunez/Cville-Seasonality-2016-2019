# ijob -p standard -A berglandlab_standard -c10 --mem=40G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(poolr)

### load test dtaa
  fl <- c("/project/berglandlab/summarized_dest_glm/glm.out.DE_Bro_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.DE_Mun_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.FI_Aka_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.PA_li_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.TR_Yes_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.UA_Ode_0.Rdata",
          "/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata"
    )

    tmp <- glm.out[mod=="aveTemp" & !is.na(rnp.clean)]


    tmp[, q:=p.adjust(p.lrt, "fdr")]


    o <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind")%do%{
      chr.tmp <- tmp[chr==chr.i]
      tab <- table(chr.tmp$q>=.1, chr.tmp$invName=="none")
      ft <- fisher.test(tab)
      data.table(chr=chr.i, est=ft$estimate, p=ft$p.value,
                invRate=mean(chr.tmp[invName!="none"]$q<=.1),
                NonInvRate=mean(chr.tmp[invName=="none"]$q<=.1), N=dim(chr.tmp)[1])

    }
    save(o, file="~/summary.Rdata")


    load("~/summary.Rdata")
    ol <- melt(o, id.vars=c("chr", "N"), measure.vars=c("invRate", "NonInvRate"))
    ol[,se:=sqrt((value*(1-value))/N)]

    ggplot(data=ol, aes(x=variable, y=value)) + geom_errorbar(aes(ymin=value-2*se, ymax=value+2*se)) + geom_point() + facet_grid(~chr)
### keep this. it works.

  o <- foreach(fli=fl)%dopar%{
    message(fli)
    load(fli)
    glm.out[mod=="aveTemp"][perm==0]
  }

  o <- rbindlist(o)
  o[,b:=as.numeric(as.character(b))]
  o[,p:=as.numeric(as.character(p))]
  o[,se:=sqrt(b^2 / qchisq(p, df=1, lower=F))]

  meta.F <- function(b.est, se){
    #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
    b.F = sum(b.est / se^2) / sum(1 / se^2)
    se.F = 1 / sqrt(sum(1 / se^2))
    p.F = pchisq( (b.F / se.F)^2, df = 1, lower = F)
    return(c(b.F, se.F, p.F, mean(sign(b.est)==1)))
  }

  o.ag <- o[!is.na(p) & !is.na(r.clean), list(est=meta.F(b.est=b, se=se), param=c("b.F", "se.F", "p.F", "st"), .N), list(chr, pos, invName)]
  o.ag <- dcast(o.ag, chr + pos + N + invName ~ param, value.var="est")
  o.ag[, r:=rank(p.F)/(length(p.F)+1)]

  fisher.test(table(o.ag[N>=1]$invName=="2Lt", o.ag[N>=1]$r<=.005))

  fisher.test(table(o.ag[N>=5]$r<=.0001, o.ag[N>=5]$st))
  table(o.ag$r<=.0001, o.ag$N)%>%prop.table()
  o.ag[!is.na(r)][order(r)]

### make composite pvalues
  o.ag <- o[!is.na(rnp.clean)][,list(cp.rnp.st=stouffer(rnp.clean)$p,
                                     cp.rnp.fi=fisher(rnp.clean)$p,
                                     cp.p.st=stouffer(p)$p,
                                     cp.p.fi=fisher(p)$p,
                                     min.rnp=min(rnp.clean), max.rnp=max(rnp.clean),
                                     min.p=min(p), max.p=max(p),
                                     mean(sign(b)),
                                     .N), list(chr, pos, invName)]
  o.ag[,rcp.rnp.st:=rank(cp.rnp.st)/(length(cp.rnp.st+1))]
  o.ag[,rcp.rnp.fi:=rank(cp.rnp.fi)/(length(cp.rnp.fi+1))]
  o.ag[,rcp.p.st:=rank(cp.p.st)/(length(cp.p.st+1))]
  o.ag[,rcp.p.fi:=rank(cp.p.fi)/(length(cp.p.fi+1))]

  o.ag[N==7][order(rcp.rnp.st)]

  o.st <- o[!is.na(rnp.clean)][,list(st=sum(sign(b[p<.05])==1), Np=sum(p<.05)), list(chr, pos)]
  setkey(o.st, chr, pos)
  setkey(o.ag, chr, pos)
  o.ag <- merge(o.ag, o.st)
  o.ag[N==7][order(rcp.rnp.st)]

### summary
  thrs <- expand.grid(sapply(c(0:9), function(x) x*10^(-5:-1)))[,1]

  o.ag.ag <- foreach(i=thrs, .combine="rbind")%dopar%{
    o.ag[,list(pr.rcp.rnp.st=mean(rcp.rnp.st<=i),
               pr.rcp.p.st=mean(rcp.p.st<=i),
               pr.rcp.rnp.fi=mean(rcp.rnp.fi<=i),
               pr.rcp.p.fi=mean(rcp.p.fi<=i),
                exp=i, .N), list(invName, chr)]

  }

  o.l <- melt(o.ag.ag, id.vars=c("invName", "chr", "exp", "N"))
  o.l[,en:=value/exp]


### windows

### define windows
  win.bp <- 10000
  step.bp <- 5000
  setkey(o.ag, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
      #chr.i <- "2L"
      tmp <- o.ag[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)

setkey(o.ag, chr, pos)
registerDoMC(5)

  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    # win.i <- 588
    i <- 0.05
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- o.ag[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
    o.tmp <- win.tmp[,list(pr=mean(r<=i), exp=i, .N)]
    o.tmp[,p:=binom.test(pr*N, N, exp)$p.value]
    o.tmp[,i:=win.i]
    return(o.tmp)
  }

win.out <- rbindlist(win.out)
win.out <- merge(win.out, wins, by="i")

### save and load
save(win.out, file="~/oL_meta_10K.Rdata")

  save(o.l, win.out, o.ag, file="~/oL.Rdata")

  library(viridis)
  load("~/oL_meta_10K.Rdata")

  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr="3R", pos=13174333/2  + 13274333/2)

### plot
  ggplot(data=o.l, aes(x=log10(exp), y=en, group=invName, color=invName)) +
  geom_line() +
  facet_grid(variable~chr)

  setnames(inv.dt, "chr.x", "chr")

  mh <- ggplot(data=win.out[N>50][order(-p)], aes(x=I(start/2 + end/2), y=-log10(p), color=pr)) +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos), color="yellow", size=2, alpha=.75) +
  geom_point() +
  facet_grid(~chr, scales="free_x") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=90))


  ggsave(mh, file="~/cp_mh_10K.png", h=4, w=10)








o.ag <- o[!is.na(rnp.clean),list(n=sum(rnp.clean<=0.05), N=length(!is.na(rnp.clean)), st=sum(sign(b)),
                z.mu=mean(z, na.rm=T)), list(chr, pos, invName)]
o.ag[,z.pr:=pnorm(z.mu, 0, 1/N)]
o.ag[N==3]
