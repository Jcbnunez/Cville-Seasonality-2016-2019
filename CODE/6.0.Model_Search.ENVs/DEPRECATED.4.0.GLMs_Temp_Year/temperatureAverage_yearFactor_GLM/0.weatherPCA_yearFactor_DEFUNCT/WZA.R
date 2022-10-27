# ijob -p standard -A berglandlab_standard -c20 --mem=20G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(tidyr)
  library(lubridate)
  library(SeqArray)
  library(metap)


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

### load GDS object
  genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### load glm output
  load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

### make windows
  glm.out <- glm.out[!is.na(rnp.clean)]
  setkey(glm.out, chr, pos)
  glm.out[,index:=c(1:dim(glm.out)[1])]
  #win.bp <- 500
  #step.bp <- 250
  #setkey(glm.out, "chr")
  #wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
  #    tmp <- glm.out[J(chr.i)]
  #    data.table(chr=chr.i,
  #                start=seq(from=min(tmp$index), to=max(tmp$index)-win.bp, by=step.bp),
  #                end=seq(from=min(tmp$index), to=max(tmp$index)-win.bp, by=step.bp) + win.bp)
  #}
  #wins[,i:=1:dim(wins)[1]]
  #dim(wins)

  ### define windows
    win.bp <- 50000
    step.bp <- 10000
    setkey(glm.out, "chr")
    wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%do%{
        tmp <- glm.out[J(chr.i)]
        data.table(chr=chr.i,
                    start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                    end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
    }
    wins[,i:=1:dim(wins)[1]]
    dim(wins)


### run windows
  setkey(glm.out, chr, pos)


  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    # win.i <- 185*5
    message(paste(win.i, dim(wins)[1], sep=" / "))
    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]
    #win.tmp <- glm.out[J(data.table(index=wins[win.i]$start:wins[win.i]$end, key="index")), nomatch=0]
    win.tmp[,Z:=qnorm(p.lrt, 0, 1)]
    win.tmp[,rnpZ:=qnorm(rnp.clean, 0, 1)]

    seqSetFilter(genofile, variant.id=unique(win.tmp$id),
              sample.id=samps[locality=="VA_ch"][year>=2014]$sampleId)

    af <- seqGetData(genofile, "annotation/format/FREQ")
    f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), id=seqGetData(genofile, "variant.id"))

    win.tmp <- merge(win.tmp, f.hat, by="id")
    win.tmp[,het:=fhat*(1-fhat)]

    win.tmp[!is.na(rnp.clean),
              list(rnp.pr=mean(rnp.clean<=0.05),
                  pr=c(0.05),
                  rbinom.p=binom.test(sum(rnp.clean<=0.05), length(rnp.clean), .05)$p.value,
                  wZa=sum(het*Z)/(sqrt(sum(het^2))),
                  wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                  rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                  rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),

                  mean.chisq=c(mean(chisq)),
                  median.chisq=c(median(chisq)),
                  gr=c("obs"),
                  nSNPs=.N, i=win.i,
                  start.bp=min(win.tmp$pos), stop.bp=max(win.tmp$pos)),
              list(mod, chr, invName, delta, locality)]
  }
  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by="i")
  win.out[mod=="aveTemp"][order(Z=wZa)]
  win.out[mod=="aveTemp"][order(rbinom.p)]


save(win.out, file="~/wza.Rdata")



### plot
  load("~/wza.Rdata")
  library(ggplot2)
  library(data.table)
  library(viridis)
  library(patchwork)
  ### load suppl data

  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)

### plot

  mh.plot.wza <- ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos), color="orange", size=2, alpha=.75) +
  geom_point(data=win.out[mod=="aveTemp"][pr==.05][order(rnp.pr)][nSNPs>50],
            aes(x=start/2 +end/2 , y=-1*log10(wZa.p),
                color=rnp.pr), size=.95) +
  facet_grid(locality~chr.x, scales="free_x") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=90))

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



  ggplot(data=win.out[mod=="aveTemp"][pr==.05][nSNPs>50],
        aes(x=1*log10(rbinom.p), y=1*log10(wZa.p))) +
  geom_point()





ggplot(win.out[mod=="aveTemp"],
  aes(x=start, y=-log10(rbinom.p), color=chr.x)) + geom_line() + facet_grid(~chr.x, scales="free_x")


ggplot(win.out[mod=="avetemp"], aes(x=-log10(wZa.p), y=log10(rbinom.p))) + geom_point()
