### libraries
  library(ggplot2)
  library(data.table)
  library(viridis)
  library(patchwork)
  library(ggridges)
  # scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest_glm_final_nested_qb/windowAnalysis/WZA_window.dt.* ~/.

### data
  load("~/WZA_window.dt.VA_ch_0.Rdata")

### load suppl data
  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr.x")

  ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)

###
  thr=0.01

### find peaks
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
     z <- i - m + 1
     z <- ifelse(z > 0, z, 1)
     w <- i + m + 1
     w <- ifelse(w < length(x), w, length(x))
     if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
   pks <- unlist(pks)
   pks
}

### 2L only
  win.chr <- win.out[mod=="aveTemp+year_factor"][pr==thr][nSNPs>100][perm==0][chr.x=="2L"]

###

  peaks.chr <- win.chr[find_peaks(x=-log10(win.chr$wZa.p), m=20)]


  mh.plot.p <-
  ggplot() +
  geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=start/1e6), size=1, alpha=.5) +
  geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=stop/1e6), size=1, alpha=.5) +
  geom_vline(data=peaks.chr, aes(xintercept=(start/2 +end/2)/1e6), color="red") +

  geom_area(data=win.chr,
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p)), size=.5, fill="lightblue") +
  geom_area(data=win.out[mod=="aveTemp+year_factor"][pr==thr][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=1*log10(rnp.binom.p)), size=.5, fill="lightgrey") +
  geom_line(data=win.chr,
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p)), size=.5, color="darkblue") +
  geom_line(data=win.out[mod=="aveTemp+year_factor"][pr==thr][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=1*log10(rnp.binom.p)), size=.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)") +
  ggtitle("Rank Normalized P\n(corrected version)")


  mh.plot.chi <-
  ggplot() +
  geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=start/1e6), size=1, alpha=.5) +
  geom_vline(data=inv.dt[chr.x=="2L"], aes(xintercept=stop/1e6), size=1, alpha=.5) +
  geom_vline(data=peaks.chr, aes(xintercept=(start/2 +end/2)/1e6), color="red") +
  geom_area(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p)), size=.5, fill="lightblue") +
  geom_area(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=1*log10(rnChi.binom.p)), size=.5, fill="lightgrey") +
  geom_line(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(wZa.p)), size=.5, color="darkblue") +
  geom_line(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
            aes(x=(start/2 +end/2)/1e6 , y=1*log10(rnChi.binom.p)), size=.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)") +
  ggtitle("Rank Normalized Chisq\n(what we've been looking at)")


  layout <- "
  AB
  "
  mega <- mh.plot.p + mh.plot.chi + plot_layout(design=layout)

  ggsave(mega , file="~/mh_plot_rnp.pdf", h=4, w=8)


#### ridgeline plot
  a1 <- ggplot(data=win.out[mod=="aveTemp+year_factor"][nSNPs>100][perm==0][chr.x=="2L"][order(-pr)]) +
  geom_ridgeline(aes(x=(start/2 +end/2)/1e6 , height=-1*log10(rnp.binom.p), y=-log10(pr)*100, group=pr), alpha=.5)


  a2 <- ggplot(data=win.out[mod=="aveTemp+year_factor"][nSNPs>100][perm==0][chr.x=="2L"][order(-pr)]) +
    geom_ridgeline(aes(x=(start/2 +end/2)/1e6 , height=log2(rnp.pr/pr), y=-log10(pr)*20, group=pr, color=as.factor(pr==0.05)), alpha=.5)

a1 + a2














  ggplot(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100], aes(x=rnp.pr, y=rnChi.pr)) +
  geom_point() +
  facet_grid(I(invName!="none")~chr.x)


    ggplot(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100],
          aes(x=-log10(rnp.binom.p), y=-log10(rnChi.binom.p))) +
    geom_point() +
    geom_abline(slope=1, intercept=0) +
    facet_grid(I(invName!="none")~chr.x)




    ggplot(data=win.out[mod=="aveTemp+year_factor"][pr==.01][order(rnp.pr)][nSNPs>100],
          aes(x=-log10(wZa.p), y=-log10(rnp.binom.p))) +
    geom_point() +
    geom_abline(slope=1, intercept=0) +
    facet_grid(I(invName!="none")~chr.x)
