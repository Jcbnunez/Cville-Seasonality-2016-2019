#system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global/crossCluster_enrichment.Rdata ~/.")

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)
  library(MASS)
  registerDoMC(4)
### data
  load("/project/berglandlab/alan/environmental_ombibus_global/crossCluster_enrichment.Rdata")
  
### get expected relative rates
  m.ag.ag.2 <- m.ag[!variable%in%c("null", "pop_year"),
                list(cross.rr.mean= mean(log2((TT/(TT+TF+FT+FF))/(thr*thr))),
                     focal.rr.mean= mean(log2((focalT/(focalT+focalF))/thr)),
                     tester.rr.mean= mean(log2((testerT/(testerT+testerF))/thr)),
                     st=mean(st.T/(st.T+st.F))),
                list(chr, inv, focalCluster, testCluster, mod, variable, thr, perm=perm)]
  m.ag.ag.2[,env:=tstrsplit(variable, "\\.")[[1]]]


#### confidence regions
### from https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in
  CIs <- foreach(chr.i=unique(m.ag.ag.2$chr), .combine="rbind")%dopar%{
    foreach(inv.i=c(T,F), .combine="rbind")%do%{
      foreach(cluster.i=unique(m.ag.ag.2$testCluster), .combine="rbind")%do%{
        foreach(thr.i=unique(m.ag.ag.2$thr), .combine="rbind")%do%{
          # chr.i <- "2L"; cluster.i<-"1.Europe_W"; inv.i<-T
          message(paste(chr.i, inv.i, cluster.i, sep=" / "))
          mv <- m.ag.ag.2[perm>0][chr==chr.i][inv==inv.i][testCluster==cluster.i][thr==thr.i]

          mv.kde <- kde2d(mv$focal.rr.mean, mv$tester.rr.mean, n = 400)
          dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
          dy <- diff(mv.kde$y[1:2])
          sz <- sort(mv.kde$z)
          c1 <- cumsum(sz) * dx * dy
          prob <- c(0.95,0.90,0.5)

          # plot:
          dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
          dc <- melt(mv.kde$z)
          dc$prob <- approx(sz,1-c1,dc$value)$y
          dc <- as.data.table(dc)
          dc[,chr:=chr.i]
          dc[,inv:=inv.i]
          dc[,testCluster:=cluster.i]
          dc[,thr:=thr.i]
          return(dc)
        }
      }
    }
  }

  setnames(CIs, c("Var1", "Var2"), c("focal.rr.mean", "tester.rr.mean"))


### plot
  p4 <-
  ggplot(data=m.ag.ag.2[perm==0][thr==.05]) +
  geom_point(aes(x=2^focal.rr.mean, y=2^tester.rr.mean, color=2^cross.rr.mean, shape=env)) +
  geom_contour(data=CIs[thr==.05], aes(x=2^focal.rr.mean, y=2^tester.rr.mean, z=prob, color=..level..), breaks=.95, color="black")+
  facet_grid(testCluster ~ chr+inv) +
  xlab("Charlottesville Relative Rate") +
  ylab("Cluster Relative Rate") +
  scale_colour_viridis(option="C", direction = -1) +
  theme_bw()

  ggsave(p4, file="~/p4_alt_3.pdf", h=7, w=14)

  m.ag.ag.2[perm==0][thr==.01][focal.rr.mean>1 & tester.rr.mean>1]
  m.ag.ag.2[perm==0][thr==.01][mod==7][chr=="2L"][inv==T][variable=="temp.max"]

  m.ag.ag.2[perm==0][thr==.01][(2^focal.rr.mean)>1 & 2^(tester.rr.mean)>1][cross.rr.mean>(2)]


#### another version

cleanMean <- function(x) median(x[x!=Inf & x!=-Inf & !is.na(x) & !is.nan(x)])
cleanSD <- function(x) sd(x[x!=Inf & x!=-Inf & !is.na(x) & !is.nan(x)])
m.ag
pad <- 0
m.ag.ag.3 <- m.ag[!variable%in%c("null", "pop_year"),
              list(cross.rr.mean=  cleanMean(log2((TT[perm==0]/(TT[perm==0]+TF[perm==0]+FT[perm==0]+FF[perm==0]) + pad) / (TT[perm!=0]/(TT[perm!=0]+TF[perm!=0]+FT[perm!=0]+FF[perm!=0]) + pad))),
                   focal.rr.mean=  cleanMean(log2((focalT[perm==0]/(focalT[perm==0]+focalF[perm==0]) + pad) / (focalT[perm!=0]/(focalT[perm!=0]+focalF[perm!=0]) + pad))),
                   tester.rr.mean= cleanMean(log2((testerT[perm==0]/(testerT[perm==0]+testerF[perm==0]) + pad) / (testerT[perm!=0]/(testerT[perm!=0]+testerF[perm!=0]) + pad))),
                    cross.rr.sd=  cleanSD(log2((TT[perm==0]/(TT[perm==0]+TF[perm==0]+FT[perm==0]+FF[perm==0]) + pad) / (TT[perm!=0]/(TT[perm!=0]+TF[perm!=0]+FT[perm!=0]+FF[perm!=0]) + pad))),
                    focal.rr.sd=  cleanSD(log2((focalT[perm==0]/(focalT[perm==0]+focalF[perm==0]) + pad) / (focalT[perm!=0]/(focalT[perm!=0]+focalF[perm!=0]) + pad))),
                   tester.rr.sd=  cleanSD(log2((testerT[perm==0]/(testerT[perm==0]+testerF[perm==0]) + pad) / (testerT[perm!=0]/(testerT[perm!=0]+testerF[perm!=0]) + pad))),
                   st.mean=cleanMean((st.pr[perm==0]/st.pr[perm!=0])),
                   st.sd=cleanSD((st.pr[perm==0]/st.pr[perm!=0]))),

              list(chr, inv, focalCluster, testCluster, mod, variable, thr)]
m.ag.ag.3[,env:=tstrsplit(variable, "\\.")[[1]]]

m.ag.ag.3[(focal.rr.mean-2*focal.rr.sd) > 0][(tester.rr.mean-2*tester.rr.sd) > 0][(cross.rr.mean-2*cross.rr.sd)>0][thr==.01]

setkey(m.ag.ag.3, mod)
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

m.ag.ag.3 <- merge(m.ag.ag.3, sets, by="mod")
m.ag.ag.3[,x:=paste(variable, mod, start, end, sep=" / ")]
tmp <- m.ag.ag.3[thr==.01][(focal.rr.mean-1.96*focal.rr.sd) > 0][(tester.rr.mean-1.96*tester.rr.sd) > 0]
table(tmp$x, tmp$inv)

[(cross.rr.mean-2*cross.rr.sd)>0]


p4 <-


ggplot(data=m.ag.ag.3[thr==.05][order(-cross.rr.mean)]) +
geom_point(data=m.ag.ag.3[thr==.05][(focal.rr.mean-2*focal.rr.sd) > 0][(tester.rr.mean-2*tester.rr.sd) > 0][(cross.rr.mean-2*cross.rr.sd)>0],
            aes(x=2^focal.rr.mean, y=2^tester.rr.mean), size=3, color="black") +
geom_point( aes(x=2^focal.rr.mean, y=2^tester.rr.mean, color=2^cross.rr.mean, shape=env)) +
facet_grid(testCluster ~ chr+inv) +
xlab("Charlottesville Relative Rate") +
ylab("Cluster Relative Rate") +
scale_colour_viridis(option="C", direction = -1) +
theme_bw()

m.ag[J(m.ag.ag.3[which.max(cross.rr.mean)])][perm==0]

m.ag.ag.3[chr=="2L"][inv==T][mod==8][variable=="temp.var"][thr==.01]

setkey(m.ag, chr, inv, focalCluster, testCluster, mod, variable, thr)
setkey(m.ag.ag.3, chr, inv, focalCluster, testCluster, mod, variable, thr)

m.ag[J(m.ag.ag.3[thr==.01][is.na(cross.rr.mean)][1])]














### old

m.ag.ag <- m.ag[,list(crossCluster.or.pr=mean(or[perm==0]>or[perm!=0], na.rm=T),
                      crossCluster.or=or[perm==0],
                      crossCluster.st.pr=mean(st.pr[perm==0]>st.pr[perm!=0], na.rm=T),
                      crossCluster.st=st.pr[perm==0],
                      focalCluster.or.pr=mean((focalT[perm==0]/focalF[perm==0]) / (focalT[perm!=0]/focalF[perm!=0]) > 1),
                      testCluster.or.pr=mean((testerT[perm==0]/testerF[perm==0]) / (testerT[perm!=0]/testerF[perm!=0]) > 1),
                      focalCluster.or=mean((focalT[perm==0]/focalF[perm==0]) / (focalT[perm!=0]/focalF[perm!=0])),
                      testCluster.or=mean((testerT[perm==0]/testerF[perm==0]) / (testerT[perm!=0]/testerF[perm!=0]))),
                list(chr, inv, focalCluster, testCluster, mod, variable, thr)]





p1 <-
ggplot(data=m.ag.ag[!variable%in%c("null", "pop_year")][thr==.01]) +
geom_tile(aes(x=variable, y=mod, fill=crossCluster.or)) +
geom_point(data=m.ag.ag[!variable%in%c("null", "pop_year")][thr==.01][crossCluster.or.pr>.99], aes(x=variable, y=mod), size=.5, color="white") +
facet_grid(chr+inv~testCluster) +coord_flip()



p2 <-
ggplot(data=m.ag.ag[!variable%in%c("null", "pop_year")][thr==.01]) +
geom_tile(aes(x=variable, y=mod, fill=crossCluster.st)) +
geom_point(data=m.ag.ag[!variable%in%c("null", "pop_year")][thr==.01][crossCluster.st.pr>.95], aes(x=variable, y=mod), size=.5, color="white") +
facet_grid(chr+inv~testCluster) +coord_flip()

mega <-
p1 / p2


ggsave(mega, file="~/crossCluster.pdf")



p3 <-
ggplot(data=m.ag.ag[!variable%in%c("null", "pop_year")][thr==.05][order(crossCluster.or.pr)]) +
geom_point(aes(x=focalCluster.or, y=testCluster.or, color=crossCluster.or.pr>.99)) +
facet_grid(testCluster ~ chr+inv)

ggsave(p3, file="~/p3.pdf")

m.ag.ag[!variable%in%c("null", "pop_year")][thr==.05][focalCluster.or>1.5 & testCluster.or>1.5][crossCluster.or.pr>.95]















m.ag.ag.2 <- m.ag[!variable%in%c("null", "pop_year"),
                list(cross.rr.mean= mean(log2((TT/(TT+TF+FT+FF))/(thr*thr))),
                     focal.rr.mean= mean(log2((focalT/(focalT+focalF))/thr)),
                     tester.rr.mean= mean(log2((testerT/(testerT+testerF))/thr)),
                     cross.rr.sd=  sd(log2((TT/(TT+TF+FT+FF))/(thr*thr))),
                     focal.rr.sd=  sd(log2((focalT/(focalT+focalF))/thr)),
                     tester.rr.sd= sd(log2((testerT/(testerT+testerF))/thr)),
                      cross.rr.uci=  quantile(log2((TT/(TT+TF+FT+FF))/(thr*thr)), .975),
                      focal.rr.uci=  quantile(log2((focalT/(focalT+focalF))/thr), .975),
                     tester.rr.uci=  quantile(log2((testerT/(testerT+testerF))/thr), .975),
                      cross.rr.lci=  quantile(log2((TT/(TT+TF+FT+FF))/(thr*thr)), .025),
                      focal.rr.lci=  quantile(log2((focalT/(focalT+focalF))/thr), .025),
                     tester.rr.lci=  quantile(log2((testerT/(testerT+testerF))/thr), .025)),
                list(chr, inv, focalCluster, testCluster, mod, variable, thr, perm=perm>0)]

m.ag.ag.3 <- m.ag.ag.2[perm==T,
                        list(focal.rr.uci=max(focal.rr.uci),
                             focal.rr.lci=min(focal.rr.lci),
                             tester.rr.uci=max(tester.rr.uci),
                             tester.rr.lci=min(tester.rr.lci)),
                        list(chr, inv, testCluster, thr)]




m.ag.ag.2[perm==0][thr==.01][focal.rr.mean>1 & tester.rr.mean>1]
