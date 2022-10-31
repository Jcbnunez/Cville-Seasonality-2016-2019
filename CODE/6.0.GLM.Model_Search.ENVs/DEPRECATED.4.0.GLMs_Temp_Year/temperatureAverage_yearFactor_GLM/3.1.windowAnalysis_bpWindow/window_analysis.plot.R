### libraries
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(viridis)
  library(patchwork)

###########################################
### load & define Inversion and Ace ROI ###
###########################################
  ### load suppl data
  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


#### load windows
  #load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/3.windowAnalysis/window_analysis_output.Rdata")
  # scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.Rdata ~/.
  #load("~/window_analysis_output.nested.Rdata")

  #scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
  load("~/window_analysis_output.nested.qb.Rdata")

### defunct
#### get sig thresholds from perms
#  win.minp.ag <- win.out[pr==0.01 & nSNPs>100,
#          list(q05=quantile(wZa.p, 0.025, na.rm=T), q5=quantile(rbinom.p, .05, na.rm=T), .N),
#          list(perm, mod, locality)]$N %>% mean
#  win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]
#
#  win.minp.ag.ag <- win.minp.ag[perm!=0, list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]
#
####


### breakpoint fix
tmp <- win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"][locality=="VA_ch"]
tmp.ag <- tmp[,.N,win.i]
sum(tmp[win.i==263]$rnp.pr * tmp[win.i==263]$nSNPs)
sum(tmp[win.i==263]$nSNPs)
rnp.pr*
### p-values
  ### summary stats v2
    win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0,
          list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
          list(mod, chr.x=chr.x, locality, win.i, start, end)]

  ### basic MH plot
    mh.plot.wza <-
    ggplot() +
    geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
    geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
    geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
    geom_ribbon(data=win.minp.ag[mod=="aveTemp+year_factor"],
                aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
                color="grey", fill="grey", alpha=0.75) +
    geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
              aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                  color=(rnp.pr)), size=.95) +
    geom_hline(yintercept=-log10(.01/1800)) +
    #geom_hline(data=win.minp.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
    facet_grid(locality~chr.x, scales="free") +
    theme_bw() +
    scale_color_viridis(option="G", direction = -1) +
    theme(axis.text.x=element_text(angle=0),
           legend.box = "vertical",
           legend.key.size=unit(1/8, 'in'),
           legend.text=element_text(size=6),
           legend.title=element_text(size=8)) +
    labs(color="Prop. top 1%", linetype="Inversion") +
    xlab("Position (Mb)") +
    ylab("-log10(Window p)")

    ggsave(mh.plot.wza, file="~/nested_qb.png", h=8.5, w=11)




### num SNPs
  ### summary stats v2

  ### basic MH plot
    mh.plot.N <-
    ggplot() +
    geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
    geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
    geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
    geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
              aes(x=(start/2 +end/2)/1e6 , y=nSNPs,
                  color=(rnp.pr)), size=.95) +
    geom_hline(yintercept=-log10(.01/1800)) +
    #geom_hline(data=win.minp.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
    facet_grid(locality~chr.x, scales="free") +
    theme_bw() +
    scale_color_viridis(option="G", direction = -1) +
    theme(axis.text.x=element_text(angle=0),
           legend.box = "vertical",
           legend.key.size=unit(1/8, 'in'),
           legend.text=element_text(size=6),
           legend.title=element_text(size=8)) +
    labs(color="Prop. top 1%", linetype="Inversion") +
    xlab("Position (Mb)") +
    ylab("nSNPs")

    ggsave(mh.plot.N, file="~/nested_qb_nSNPs.png", h=8.5, w=11)








### density plots using ggplot geom_density
  pr.densPlot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=100][pr==0.01],
            aes(x=qlogis(rnp.pr),
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=100][pr==0.01],
              aes(x=qlogis(rnp.pr),
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod+as.factor(invName=="none")) +
  theme_bw() +
  theme(legend.position="none")


  wza.densPlot <-
  ggplot() +
  #geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=100][pr==0.01],
            aes(x=-wZa,
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=100][pr==0.01],
              aes(x=-wZa,
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod+as.factor(invName=="none")) +
  theme_bw() +
  theme(legend.position="none")


  mega <-
  pr.densPlot + wza.densPlot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="~/densPlot_nested.qb.split.pdf", h=10, w=10)


###  quantile plot

  win.out.quan <- win.out[nSNPs>100,
                          list(rnp.lci=quantile(rnp.pr, 0.025), rnp.uci=quantile(rnp.pr, 0.975), rnp.med=quantile(rnp.pr, 0.5)),
                          list(chr.x, inv=as.factor(invName!="none"), mod, perm, pr, locality)]

  win.out.quan.ag <- win.out.quan[,
                              list(rnp.lci=mean(qlogis(rnp.lci)), rnp.uci=mean(qlogis(rnp.uci)), rnp.med=mean(qlogis(rnp.med))),
                              list(chr.x, inv, mod, gr=(perm==0), pr, locality)]



  ggplot(data=win.out.quan.ag[pr==0.01][locality=="VA_ch"]) +
  geom_point(aes(y=inv, x=rnp.med, color=gr)) +
  geom_linerange(aes(y=inv, xmin=rnp.lci, xmax=rnp.uci, color=gr)) +

  facet_grid(mod~chr.x)




### noodling around
  ggplot(data=win.out[nSNPs>50][locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"][perm<10]) +
  geom_point(aes(x=qlogis(rnp.pr), y=-log10(rbinom.p), size=nSNPs)) +
  facet_wrap(~perm)

  ggplot() +
  geom_point(data=win.out[nSNPs>50][locality=="VA_ch"][mod=="aveTemp+year_factor"][order(perm)][perm!=0],
            aes(x=rnp.pr, y=-log10(rbinom.p), group=perm), color="grey", alpha=.5) +
  geom_point(data=win.out[nSNPs>50][locality=="VA_ch"][mod=="aveTemp+year_factor"][order(perm)][perm==0],
            aes(x=rnp.pr, y=-log10(rbinom.p), group=perm), color="red") +
  facet_grid(as.factor(invName!="none")~chr.x)




library(patchwork)
### more noodling(
wza <- ggplot(data=win.out[mod=="aveTemp+year_factor"][locality=="VA_ch"][pr==0.01][perm==0][chr.x=="2L"],
        aes(x=I(start/2 + end/2), y=-log10(rnp.wZa.p))) +

               geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=start, linetype=invName)) +
               geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=stop, linetype=invName)) +
geom_vline(aes(xintercept=c(9891255/2 + 9893446/2)), color="red", size=2) +
geom_vline(aes(xintercept=c(5100877/2 + 5207002/2)), color="green", size=2) +
geom_hline(yintercept=-log10(.01/1800)) +
geom_line()

rbin <- ggplot(data=win.out[mod=="aveTemp+year_factor"][locality=="VA_ch"][pr==0.01][perm==0][chr.x=="2L"],
        aes(x=I(start/2 + end/2), y=-log10(rbinom.p))) +

               geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=start, linetype=invName)) +
               geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=stop, linetype=invName)) +
geom_vline(aes(xintercept=c(9891255/2 + 9893446/2)), color="red", size=2) +
geom_vline(aes(xintercept=c(5100877/2 + 5207002/2)), color="green", size=2) +
geom_hline(yintercept=-log10(.01/1800)) +
geom_line()



win.out[mod=="aveTemp+year_factor"][locality=="VA_ch"][pr==0.01][perm==0][chr.x=="2L"][wZa.p<1e-120]

p1 <- wza / rbin

ggsave(p1, file="~/msp_mthelicase_lineup.pdf")

       )




       geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
       geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
       geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
       geom_ribbon(data=win.minp.ag[mod=="aveTemp+year_factor"],
                   aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
                   color="grey", fill="grey", alpha=0.75) +
       geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
                 aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                     color=(rnp.pr)), size=.95) +
       geom_hline(yintercept=-log10(.01/1800)) +




,002 [+]

5100877/2 + 5207002/2
