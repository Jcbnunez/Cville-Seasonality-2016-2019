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
  load("~/window_analysis_output.nested.Rdata")

### get sig thresholds from perms

  win.minp.ag <- win.out[wZa.p>0 & nSNPs>50,
                        list(minp.wza=min(wZa.p, na.rm=T), minp.rnp=min(rbinom.p, na.rm=T)),
                        list(perm, mod, chr.x=chr.x, locality)]
  win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp"]

  win.minp.ag.ag <- win.minp.ag[, list(q5.minp.wza=quantile(minp.wza[perm!=0], .05, na.rm=T),
                                       q5.minp.rnp=quantile(minp.rnp[perm!=0], .05, na.rm=T),
                                       obs.minp.wza=minp.wza[perm==0],
                                       obs.minp.rnp=minp.rnp[perm==0]),
                                  list(locality, mod, chr.x)]

  win.minp.ag.ag[obs.minp.rnp<q5.minp.rnp]
  win.minp.ag.ag[obs.minp.wza<q5.minp.wza]

### basic MH plot
  mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>50][perm==0],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p),
                color=rnp.pr), size=.95) +
  geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(q5.minp.wza))) +
  facet_grid(locality~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
         legend.box = "vertical",
         legend.key.size=unit(1/8, 'in'),
         legend.text=element_text(size=6),
         legend.title=element_text(size=8)) +
  labs(color="Prop. top 5%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)")

  ggsave(mh.plot.wza, file="~/mh_plot_wza_nested.pdf")




### distribution plots

  foreach(pop=unique(win.out$locality))%do%{
    #pop <- "UA_Ode"
    message(pop)

    pr.densPlot <-
    ggplot() +
    geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
    geom_density(data=win.out[locality==pop][perm!=0][nSNPs>=50],
              aes(x=qlogis(rnp.pr),
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="grey") +
    geom_density(data=win.out[locality==pop][perm==0][nSNPs>=50],
                aes(x=qlogis(rnp.pr),
                    group=interaction(perm, (invName=="none")),
                    linetype=as.factor(invName=="none")),
                color="black", size=1) +
    facet_grid(chr.x~mod) +
    theme_bw() + theme(axis.text=element_text(size=10)) +
    theme(legend.position="bottom") +
    xlab("Proporiton of SNPs in Window with RNP < 0.01; logit scale")


    wza.densPlot <-
    ggplot() +
    #geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
    geom_density(data=win.out[locality==pop][perm!=0][nSNPs>=50],
              aes(x=-wZa,
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="grey") +
    geom_density(data=win.out[locality==pop][perm==0][nSNPs>=50],
                aes(x=-wZa,
                    group=interaction(perm, (invName=="none")),
                    linetype=as.factor(invName=="none")),
                color="black", size=1) +
    facet_grid(chr.x~mod) +
    theme_bw() + theme(axis.text=element_text(size=10)) +
    theme(legend.position="bottom") +
    xlab("Window Z statistic")


    mega <-
    pr.densPlot + wza.densPlot +
    plot_annotation(tag_levels="A", title=pop)

    ggsave(mega, file=paste("~/densPlot_nested.", pop, ".png", sep=""), h=10, w=10)

}



summary(lm(rnp.pr~nSNPs, win.out[locality=="VA_ch"][perm!=0]))




  std <- density(qlogis(win.out[locality=="VA_ch"][perm==0][nSNPs>=50][chr.x=="2L"][invName=="none"][mod=="aveTemp+year_factor"]$rnp.pr), from=-8, to=8)
  inv <- density(qlogis(win.out[locality=="VA_ch"][perm==0][nSNPs>=50][chr.x=="2L"][invName=="2Lt"][mod=="aveTemp+year_factor"]$rnp.pr) , from=-8, to=8)

  ggplot() +
  geom_line(aes(x=std$x, y=std$y), color="red") +
  geom_line(aes(x=inv$x, y=inv$y), color="black")




    ggsave(densPlot, file="~/densPlot.pdf")
