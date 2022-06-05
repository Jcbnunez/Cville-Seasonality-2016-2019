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


### distribution plots
  dFun <- function(x, min, max) {
    dens <- density(x, from=min, to=max)
    c(paste("x_", dens$x, sep=""), paste("y_", dens$y, sep=""))
  }


### RNP distributions
  dens.ag <- win.out[,
                      list(d=dFun(qlogis(rnp.pr), min=-6, max=6)),
                      list(chr.x=chr.x, inv=invName!="none", mod, locality, perm)]

  dens.ag[,value:=tstrsplit(d, "_")[[2]]%>%as.numeric]
  dens.ag[,variable:=tstrsplit(d, "_")[[1]]]


  dens.ag.x <- dens.ag[variable=="x"][,-c("d", "variable"),with=F]
  dens.ag.y <- dens.ag[variable=="y"][,-c("d", "variable"),with=F]

  rnp.dens.ag.w <- dens.ag.x
  setnames(rnp.dens.ag.w, "value", "x")
  rnp.dens.ag.w[,y:=dens.ag.y$value]
  rnp.dens.ag.w[,stat:="rnp"]

### wZa distributions
  dens.ag <- win.out[wZa!=Inf & wZa!=-Inf,
                      list(d=dFun(wZa, min=-60, max=10)),
                      list(chr.x=chr.x, inv=invName!="none", mod, locality, perm)]

  dens.ag[,value:=tstrsplit(d, "_")[[2]]%>%as.numeric]
  dens.ag[,variable:=tstrsplit(d, "_")[[1]]]


  dens.ag.x <- dens.ag[variable=="x"][,-c("d", "variable"),with=F]
  dens.ag.y <- dens.ag[variable=="y"][,-c("d", "variable"),with=F]

  wza.dens.ag.w <- dens.ag.x
  setnames(wza.dens.ag.w, "value", "x")
  wza.dens.ag.w[,y:=dens.ag.y$value]
  wza.dens.ag.w[,stat:="wZa"]


### plot components
  pr.densPlot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_line(data=rnp.dens.ag.w[locality=="VA_ch"][perm!=0],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="grey") +
  geom_line(data=rnp.dens.ag.w[locality=="VA_ch"][perm==0],
              aes(x=x, y=y,
                  group=interaction(perm, inv),
                  linetype=inv),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")


  wza.densPlot <-
  ggplot() +
  #geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_line(data=wza.dens.ag.w[locality=="VA_ch"][perm!=0],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="grey") +
  geom_line(data=wza.dens.ag.w[locality=="VA_ch"][perm==0],
              aes(x=x, y=y,
                  group=interaction(perm, inv),
                  linetype=inv),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")


  mega <-
  pr.densPlot + wza.densPlot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="~/densPlot_nested_v2.pdf", h=10, w=10)




#### version that uses counts
  win.out[,rnp.pr.bin:=round(qlogis(rnp.pr),1)]
  win.out[,wZa.bin:=round(wZa*2)/2]

  rnp.ag <- win.out[,list(y=.N), list(mod, chr.x, inv=invName!="none", locality, perm, x=rnp.pr.bin)]
  wza.ag <- win.out[,list(y=.N), list(mod, chr.x, inv=invName!="none", locality, perm, x=wZa.bin)]

  rnp.plot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_line(data=rnp.ag[locality=="VA_ch"][perm!=0],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="grey") +
  geom_line(data=rnp.ag[locality=="VA_ch"][perm==0],
              aes(x=x, y=y,
                  group=interaction(perm, inv),
                  linetype=inv),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")


  wza.plot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_line(data=wza.ag[locality=="VA_ch"][perm!=0],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="grey") +
  geom_line(data=wza.ag[locality=="VA_ch"][perm==0],
              aes(x=x, y=y,
                  group=interaction(perm, inv),
                  linetype=inv),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")



  mega <-
  pr.densPlot + wza.densPlot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="~/densPlot_nested_v3.pdf", h=10, w=10)




























  pr.densPlot <-
  ggplot() +
  geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=50],
            aes(x=qlogis(rnp.pr),
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=50],
              aes(x=qlogis(rnp.pr),
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")


  wza.densPlot <-
  ggplot() +
  #geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
  geom_density(data=win.out[locality=="VA_ch"][perm!=0][nSNPs>=50],
            aes(x=-wZa,
                group=interaction(perm, (invName=="none")),
                linetype=as.factor(invName=="none")),
            color="grey") +
  geom_density(data=win.out[locality=="VA_ch"][perm==0][nSNPs>=50],
              aes(x=-wZa,
                  group=interaction(perm, (invName=="none")),
                  linetype=as.factor(invName=="none")),
              color="black", size=1) +
  facet_grid(chr.x~mod) +
  theme_bw() +
  theme(legend.position="none")


  mega <-
  pr.densPlot + wza.densPlot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="~/densPlot_nested.pdf", h=10, w=10)





summary(lm(rnp.pr~nSNPs, win.out[locality=="VA_ch"][perm!=0]))




  std <- density(qlogis(win.out[locality=="VA_ch"][perm==0][nSNPs>=50][chr.x=="2L"][invName=="none"][mod=="aveTemp+year_factor"]$rnp.pr), from=-8, to=8)
  inv <- density(qlogis(win.out[locality=="VA_ch"][perm==0][nSNPs>=50][chr.x=="2L"][invName=="2Lt"][mod=="aveTemp+year_factor"]$rnp.pr) , from=-8, to=8)

  ggplot() +
  geom_line(aes(x=std$x, y=std$y), color="red") +
  geom_line(aes(x=inv$x, y=inv$y), color="black")




    ggsave(densPlot, file="~/densPlot.pdf")
