### libraries
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(viridis)
  library(patchwork)


#### load windows
  #load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/3.windowAnalysis/window_analysis_output.Rdata")
  # scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.Rdata ~/.
  load("~/window_analysis_output.nested.Rdata")


### summarize
  win.out.orig <- win.out[perm==0]
  win.out.perm <- win.out[perm!=0]

  setkey(win.out.orig, i, mod, chr.x, invName, locality, delta)
  setkey(win.out.perm, i, mod, chr.x, invName, locality, delta)

  m <- merge(win.out.orig, win.out.perm)

  win.out.ag <- m[nSNPs.x>50 & nSNPs.y>50 & rnp.pr.x>0 & rnp.pr.x<1 & rnp.pr.y>0 & rnp.pr.y<1,
                        list(t=t.test(qlogis(rnp.pr.x), qlogis(rnp.pr.y), paired=T)$statistic),
                        list(mod, chr.x, inv=as.factor(invName!="none"), locality, perm=perm.y)]


  ggplot(data=win.out.ag) +
  geom_boxplot(aes(x=chr.x, y=t, fill=inv, group=interaction(inv, chr.x)), position="dodge") +
  facet_grid(mod~locality)


interaction(chr.x, invName!="none")
