scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest_glm_morePerms_nested_qb/windowAnalysis_lotsPR_snpWin/WZA_window.dt.VA_ch_0.Rdata ~/.


### libraries
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(viridis)
  library(patchwork)

### load
  load("/Users/alanbergland/WZA_window.dt.VA_ch_0.Rdata")


  ###########################################
  ### load & define Inversion and Ace ROI ###
  ###########################################
    ### load suppl data
    inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
    setnames(inv.dt, "chrom", "chr.x")

    ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


  ### basic MH plot
    mh.plot.wza <-
    ggplot() +
    geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
    geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
    geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
    geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
              aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                  color=(rnp.pr)), size=.95) +
    geom_hline(yintercept=-log10(.05/2845)) +
    #geom_hline(data=win.minp.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
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

    ggsave(mh.plot.wza, file="~/nested_qb.png", h=8.5, w=11)


### summary
  max(win.out$i)

  win.out.ag <- win.out[,list(rr=mean(log2(rnp.pr/pr)),
                              rr.sd=sd(log2(rnp.pr/pr))/sqrt(length(rnp.pr))),
                        list(pr, mod, chr=chr.x, inv=invName!="none")]

  ggplot(data=win.out.ag) +
  geom_hline(yintercept=1) +
  geom_ribbon(aes(x=log10(pr), ymin=(rr-2*rr.sd), ymax=(rr+2*rr.sd), group=mod, fill=mod), alpha=.5)+
  geom_line(aes(x=log10(pr), y=(rr), group=mod, color=mod)) +
  facet_grid(chr~inv)



  ggplot(data=win.out[mod=="aveTemp+year_factor"][pr==.01]) +
  geom_point(aes(x=snpDensity, y=-log10(rbinom.p)))

  win.out[chr.x=="2L"][start>=9450000][end<=10000000][pr==0.05][mod=="aveTemp+year_factor"]
  win.out[i>200 & i<210][pr==0.05][mod=="aveTemp+year_factor"]



  win.out[,snpDensity:=nSNPs/abs(start-end)]
