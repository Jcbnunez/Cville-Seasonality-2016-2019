# ijob -p standard -A berglandlab_standard -c20 --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  #registerDoMC(20)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(SeqArray)
  library(lubridate)

### run windows
  fl <- list.files("/scratch/aob2x/summarized_dest_glm/", "WZA_window", full.names=T)
  fl <- fl[!grepl("qb", fl)]
  o <- foreach(fl.i=fl, .combine="rbind")%do%{
    #fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    #win.out[,psi:=tstrsplit(fl.i, "/") %>% last ]
    #win.out[,psi:=tstrsplit(psi, "\\.")[[2]]]
    win.out[,perm:=tstrsplit(psi, "_")[[3]]]
    win.out[,locality:=paste(tstrsplit(psi, "_")[[1]], tstrsplit(psi, "_")[[2]], sep="_")]

    win.out
  }


  save(o, file="~/Overwintering_18_19/GEA_megaplot/windows_wza.Rdata")


### plot
  library(ggplot2)
  library(data.table)
  library(viridis)
  library(patchwork)

### scp aob2x@rivanna.hpc.virginia.edu:~/windows_wza.Rdata ~/.
  load("~/windows_wza.qb.Rdata")

### summarize the perms
  o.minp.ag <- o[wZa.p>0, list(minp=min(wZa.p, na.rm=T)), list(perm, mod, chr.x=chr.x, locality)]
  o.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp"]
  o.minp.ag[locality=="DE_Mun"][chr.x=="2L"][mod=="aveTemp"]

  o.minp.ag.ag <- o.minp.ag[perm!=0, list(q5.minp=quantile(minp, .025, na.rm=T)), list(locality, mod, chr.x)]
  o.minp.ag.ag[mod=="aveTemp"][locality=="DE_Mun"]

### reload and subset to Ace to ground truth

#  o.minp.ag <- o[, list(minp=min(wZa.p, na.rm=T)), list(perm, mod, chr=chr.x, locality)]
#  o.minp.ag[locality=="VA_ch"][chr=="2L"][mod=="aveTemp"]
#  o.minp.ag[locality=="FI_Aka"][chr=="2L"][mod=="aveTemp"]
#
#  o.minp.ag.ag <- o.minp.ag[perm!=0, list(min.minp=min(minp, na.rm=T)), list(locality, mod, chr)]
#  o.minp.ag.ag
#


### load suppl data
  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)

### plot ##[locality%in%c("TR_Yes", "UA_Ode", "VA_ch")]
  mh.plot.wza <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_point(data=o[mod=="aveTemp"][pr==.05][order(rnp.pr)][nSNPs>50][perm==0],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p),
                color=rnp.pr), size=.95) +
  geom_hline(data=o.minp.ag.ag[mod=="aveTemp"], aes(yintercept=-log10(q5.minp))) +
  facet_grid(locality~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0)) +
  labs(color="Prop. top 5%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)")

  ggsave(mh.plot.wza, file="~/mh.plot.wza.qb.png", h=8, w=11)

### window meta-analysis plot #[locality%in%c("TR_Yes", "UA_Ode", "VA_ch")]
  source("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/metap.R")
  om <- o[wZa!=Inf & wZa!=-Inf,
                   list(win.metaZ=sum(wZa)/sqrt(.N),
                        win.metaZ.p=pnorm(sum(wZa)/sqrt(.N)),
                        win.rnp.pr=sum(rnp.pr*nSNPs)/sum(nSNPs),
                        nLocales=.N,
                        nSNPs=mean(nSNPs)) ,
                   list(perm, win.i=i, chr.x=chr.x, start, end, mod, invName)]
  om[mod=="aveTemp"][order(win.metaZ)]
  om[mod=="aveTemp"][order(win.rnp.pr)]

  om.ag <- om[perm!=0,list(lci=quantile(win.metaZ.p, .025),
                           uci=quantile(win.metaZ.p, .975)),
                      list(chr.x, mod)]


  meta.plot <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_point(data=om[mod=="aveTemp"][perm==0][nSNPs>50],
              aes(x=(start/2 + end/2)/1e6, y=-log10(win.metaZ.p),
                  color=win.rnp.pr)) +
  geom_hline(data=om.ag[mod=="aveTemp"], aes(yintercept=-log10(lci))) +
  facet_grid(~chr.x, scales="free_x") +
  scale_color_viridis(option="G", direction = -1) +
  theme_bw() +
  labs(color="Prop. top 5%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window meta-p)")

  ggsave(meta.plot, file="~/meta_plot.qb.png")


### cville only
  mh.plot.wza.cville <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_point(data=o.small[mod=="aveTemp"][pr==.05][order(rnp.pr)][nSNPs>50][locality=="VA_ch"],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p),
                color=rnp.pr), size=1.5) +
  geom_hline(data=o.minp.ag.ag[mod=="aveTemp"][locality=="VA_ch"], aes(yintercept=-log10(q5.minp))) +
  facet_grid(locality~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0)) +
  labs(color="Prop. top 5%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)")

  ggsave(mh.plot.wza.cville, file="~/mh.plot.wza.cville.png")






  o[perm==9][start==30074333 & end==30174333]











  mh.plot.binom <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos), color="orange", size=2, alpha=.75) +
  geom_point(data=win.out[mod=="aveTemp"][pr==.05][order(rnp.pr)][nSNPs>50],
            aes(x=start/2 +end/2 , y=-1*log10(rnp.wZa.p),
                color=rnp.pr), size=.95) +
  facet_grid(locality~chr.x, scales="free_x") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=90))


  #+ ylim(-25, 150) +

  mega <- mh.plot.wza / mh.plot.binom
  ggsave(mega, file="~/mh_plot_wza.png", h=8, w=10.5)

























#o2 <- o[,list(rnp.pr=sum(rnp.pr*nSNPs)/sum(nSNPs), nSNPs=sum(nSNPs)), list(i, mod, chr.x, delta, gr, start, end, psi, perm, pr, locality)]


#o.ag <- o2[,list(rnpr.pr=rnp.pr[perm==0], nSNPs=nSNPs[perm==0],
#                lci.pr=quantile(rnp.pr[perm!=0], .1),  uci.pr=quantile(rnp.pr[perm!=0], .9),
#                p.binom=binom.test(rnp.pr[perm==0] * nSNPs[perm==0], nSNPs[perm==0], pr)$p.value),
#          list(chr.x, start, end, mod, pr, locality)]

### load in allele frequencies and metadata
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


### open GDS for common SNPs (PoolSNP)
 genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### extract allele frequencies function
  getData <- function(v.id, samplesUse=samps$sampleId) {
    # v.id=glm.out[chr=="3R"][pos>=13174333][pos<=13274333][rnp.clean<.001][mod=="aveTemp"]$id
    # samplesUse=samps[locality=="VA_ch"][year>2012]$sampleId

    ### filter to target
      seqSetFilter(genofile, variant.id=v.id, sample.id=samplesUse, verbose=T)

    ### get frequencies
      message("Allele Freqs")

      ad <- seqGetData(genofile, "annotation/format/AD")
      dp <- seqGetData(genofile, "annotation/format/DP")

      #if(class(dp)[1]!="SeqVarDataList") {
      #  dp.list <- list()
      #  dp.list$data <- dp
      #  dp <- dp.list
      #}

      af <- data.table(chr=rep(seqGetData(genofile, "chromosome"), each=dim(ad$data)[1]),
                       pos=rep(seqGetData(genofile, "position"), each=dim(ad$data)[1]),
                       ad=expand.grid(ad$data)[,1],
                       dp=expand.grid(dp$data)[,1],
                       sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                       variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
      message("merge")
      #afi <- merge(af, snp.dt1.an, by="variant.id")
      #afi <- merge(af, snp.dt, by.x="variant.id", by.y="id")
      afi <- af
      afi[,af:=ad/dp]

    ### calculate effective read-depth
      afis <- merge(afi, samps, by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
      afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      #afis[,-c("n"), with=F]
      #afis[,c("sampleId", "af_nEff", "nEff"), with=F]
      afis
  }


  save(afis, file="~/afis.Rdata")

  load("~/afis.Rdata")
  library(ggplot2)
  library(data.table)
  ggplot(data=afis, aes(x=aveTemp, y=af_nEff, group=variant.id)) + geom_line() + geom_smooth(method = "glm", se = FALSE)

  summary(glm(af_nEff~aveTemp+as.factor(variant.id), data=afis, weights=afis$nEff, family=binomial())

  o.ag[locality=="VA_ch"][p.binom<1e-30][mod=="aveTemp"]

  load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

  glm.out[chr=="3R"][pos>=13174333][pos<=13274333][rnp.clean<.05][mod=="aveTemp"]








### scp aob2x@rivanna.hpc.virginia.edu:~/windows.Rdata ~/.

library(data.table)
library(ggplot2)

load("~/windows.Rdata")
o2 <- o[,list(rnp.pr=sum(rnp.pr*nSNPs)/sum(nSNPs), nSNPs=sum(nSNPs)), list(i, mod, chr.x, delta, gr, start, end, psi, perm, pr, locality)]

inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
setnames(inv.dt, "chrom", "chr.x")

ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


ggplot() +
geom_vline(data=inv.dt, aes(xintercept=start, color=invName)) +
geom_vline(data=inv.dt, aes(xintercept=stop, color=invName)) +
geom_line(data=o[mod=="aveTemp"][perm==0][pr==.05][nSNPs>50], aes(x=start/2 +end/2 , y=rnp.pr)) +
facet_grid(locality~chr.x)





o.ag <- o2[,list(rnpr.pr=rnp.pr[perm==0], nSNPs=nSNPs[perm==0],
                lci.pr=quantile(rnp.pr[perm!=0], .1),  uci.pr=quantile(rnp.pr[perm!=0], .9),
                p.binom=binom.test(rnp.pr[perm==0] * nSNPs[perm==0], nSNPs[perm==0], pr)$p.value),
          list(chr.x, start, end, mod, pr, locality)]


mh.plot <- ggplot() +
geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
geom_vline(data=ace.dt, aes(xintercept=pos), color="orange", size=2, alpha=.75) +
geom_point(data=o.ag[mod=="aveTemp"][pr==.05][nSNPs>50][order(-rnpr.pr)],
          aes(x=start/2 +end/2 , y=-1*sign(rnpr.pr - pr)*log10((p.binom)),
              color=rnpr.pr), size=.95) +
facet_grid(locality~chr.x, scales="free_x") +
theme_bw() +
scale_color_viridis(option="G", direction = -1) +
theme(axis.text.x=element_text(angle=90))

#+ ylim(-25, 150) +
ggsave(mh.plot, file="~/mh_plot_temp.png", h=4, w=10.5)







mh.plot <- ggplot() +
geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
geom_vline(data=ace.dt, aes(xintercept=pos), color="orange", size=2, alpha=.95) +
geom_point(data=o.ag[mod=="aveTemp"][pr==.05][nSNPs>250][order(-rnpr.pr)],
          aes(x=start/2 +end/2 , y=rnpr.pr,
              color=-1*log10(p.adjust(p.binom))), size=.95) +
geom_ribbon(data=o.ag[mod=="aveTemp"][pr==.05][nSNPs>50],
            aes(x=start/2 +end/2, ymin=uci.pr, ymax=lci.pr), fill="black", alpha=.5) +

facet_grid(locality~chr.x, scales="free_x") +
theme_bw() +
scale_color_viridis(option="A", direction = -1) +
theme(axis.text.x=element_text(angle=90))

#+ ylim(-25, 150) +
ggsave(mh.plot, file="~/mh_plot.png", h=8, w=10.5)

o.ag[nSNPs>250][chr.x=="2R"][pr==0.05][mod=="aveTemp"][order(rnpr.pr)][locality=="VA_ch"]








o.ag.ag <- o.ag[mod=="aveTemp"][pr==.05][nSNPs>250][,list(n=sum(rnpr.pr>.10)),
          list(chr.x, start, end, mod, pr)]

o.ag.ag[n>2]

count.plot <- ggplot() +
geom_vline(data=inv.dt, aes(xintercept=start, color=invName)) +
geom_vline(data=inv.dt, aes(xintercept=stop, color=invName)) +
geom_vline(data=ace.dt, aes(xintercept=pos), color="yellow", size=2, alpha=.5) +
geom_line(data=o.ag.ag[mod=="aveTemp"][pr==.05], aes(x=start/2 +end/2 , y=n)) +
facet_grid(.~chr.x, scales="free_x") + ylim(-0, 25) +
geom_vline(data=ace.dt, aes(xintercept=pos), color="yellow", size=2, alpha=.5) +
theme_bw()
ggsave(count.plot, file="~/mh_count.png", h=8, w=10.5)


dim(glm.out[mod=="aveTemp"][perm==0][chr=="2L"][pos==6188948])




fl <- c("/project/berglandlab/summarized_dest_glm/glm.out.UA_Ode_0.Rdata",
        "/project/berglandlab/summarized_dest_glm/glm.out.DE_Mun_0.Rdata",
         "/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

o <- foreach(fli=fl)%dopar%{
  message(fli)
  load(fli)
  glm.out[mod=="aveTemp"][perm==0]
}

o <- rbindlist(o)
o[,b:=as.numeric(as.character(b))]
o[,z:=qnorm(rnp.clean, 0, 1)]

o.ag <- o[!is.na(rnp.clean),list(n=sum(rnp.clean<=0.05), N=length(!is.na(rnp.clean)), st=sum(sign(b)),
                z.mu=mean(z, na.rm=T)), list(chr, pos, invName)]
o.ag[,z.pr:=pnorm(z.mu, 0, 1/N)]
o.ag[N==3]

### define windows
  win.bp <- 50000
  step.bp <- 10000
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
registerDoMC(10)
win.out <- foreach(win.i=c(1:dim(wins)[1]))%dopar%{
  # win.i <- 954
  message(paste(win.i, dim(wins)[1], sep=" / "))
  win.tmp <- o.ag[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
  o.tmp <- win.tmp[N==3][,list(pr=sum(n==3, na.rm=T), meanZ=mean(z.mu, na.rm=T),
                               N=length(n[!is.na(n)]),
                               exp=(.05)^3*length(n[!is.na(n)]),
                               mean.beta.st=mean(st, na.rm=T),
                               n3.beta.st=mean(st[n==3], na.rm=T))]
  o.tmp[,i:=win.i]
  return(o.tmp)
}

win.out <- rbindlist(win.out)

win.out <- merge(win.out, wins, by="i")
save(win.out, o.ag, o, file="~/3pop_win.Rdata")

scp aob2x@rivanna.hpc.virginia.edu:~/3pop_win.Rdata ~/.

library(data.table)
library(ggplot2)
load("~/3pop_win.Rdata")

ggplot(data=win.out[N>50], aes(x=i, y=log2((1+pr)/(1+exp)))) + geom_line()
ggplot(data=win.out[N>50], aes(x=i, y=n3.beta.st, color=chr)) + geom_line()
ggplot(data=win.out[N>50], aes(x=i, y=meanZ, color=chr)) + geom_line()

ggplot(data=win.out[N>50], aes(x=log2((1+pr)/(1+exp)), y=mean.beta.st, color=chr)) + geom_point()
summary(lm(mean.beta.st~log2((1+pr)/(1+exp))*chr, win.out[N>50]))

ggplot(data=o.ag[N==3][z.pr<.05], aes(x=invName, y=-log10(z.pr), group=invName, color=invName)) + geom_boxplot() + facet_grid(~chr)

ggplot(data=o.ag[N==3][z.pr<.05], aes(x=pos, y=-log10(z.pr), group=chr, color=chr)) + geom_point() + facet_grid(~chr)


o.ag[,z.pr.r := rank(z.pr)/(length(z.pr)+1)]
o.ag[,list(pr=mean(z.pr.r<.001)), list(chr, invName)]


win.out[N>50][(pr/exp)>275]

o.ag[chr=="2L"][N==3][,list(pr=sum(n==3, na.rm=T), .N, exp=(.05)^3*.N), list(smal=I(pos>=6179133 & pos<=6211114))]


o.ag[N==3][,list(pr=sum(n==3, na.rm=T), .N, exp=(.05)^3*.N), list(invName)]




dim(glm.out[mod=="aveTemp"][perm==0][chr=="2L"][pos>=6179133 & pos<=6211114])
360*.05


o.ag[start<=13154180 & end>=13154180][chr.x=="2L"][pr==0.05][mod=="aveTemp"]
o[start<=13154180 & end>=13154180][chr.x=="2L"][pr==0.05][mod=="aveTemp"][locality=="VA_ch"][perm==0]
o[end>=13154180][chr.x=="2L"][pr==0.05][mod=="aveTemp"][locality=="VA_ch"][i=262]
o.ag[chr.x=="3R"][pr==0.05][mod=="aveTemp"][locality=="VA_ch"][p.binom<1e-25]


o.ag1 <- o[,list( rnpr.pr=rnp.pr,
                p.binom=binom.test(rnp.pr * nSNPs, nSNPs, pr)$p.value),
          list(chr.x, start, end, mod, pr, locality, perm)]



  mh.plot <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_line(data=o.ag1[mod=="aveTemp"][pr==.05][perm==2 | perm==0],
          aes(x=start/2 +end/2 ,
              y=-1*sign(rnpr.pr - pr)*log10(p.adjust(p.binom)),
              color=as.factor(perm==0))) +
  facet_grid(locality~chr.x) + ylim(-100, 100)

  ggsave(mh.plot, file="~/mh_plot.png", h=8, w=10.5)






o.ag[mod=="aveTemp"][chr.x=="2L"][start==13055437][pr==0.05]

ggplot() +
geom_vline(data=inv.dt, aes(xintercept=start, color=invName)) +
geom_vline(data=inv.dt, aes(xintercept=stop, color=invName)) +
geom_line(data=o.ag[mod=="aveTemp"][pr==.05], aes(x=start/2 +end/2 , y=(rr))) +
facet_grid(locality~chr.x) + ylim(.9, 1.2)













  o.ag2.ag <- o.ag2[,list(or=mean((or)),
                  or2=mean((or2))),

              list(model, factor, locality)]



ggplot(data=o.ag2.ag) +
 geom_point(aes(x=locality, y=or)) +
 geom_point(aes(x=locality, y=or2), color="red") +
 facet_grid(model~factor) +
 theme(axis.text.x=element_text(angle=90))
