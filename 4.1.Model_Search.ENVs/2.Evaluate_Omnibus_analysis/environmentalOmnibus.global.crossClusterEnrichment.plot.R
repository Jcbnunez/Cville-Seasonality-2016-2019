#system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global/crossCluster_enrichment.Rdata ~/.")

### libraries
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(MASS)
  library(foreach)
  library(doMC)
  library(viridis)
  registerDoMC(4)
  library(doParallel)
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(magrittr)
  library(doMC)
  library(SeqArray)
  library(lubridate)

  sets <- data.table(mod=c(1:11),
                     start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                     end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
  sets
  
  
  ### open GDS
  genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
  
  ### samps
  samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps <- samps[set!="dgn"]
  samps[,Date:=date(paste(year, month, day, sep="-"))]
  
  ###
  base <- "/project/berglandlab/alan/environmental_ombibus_global"
  
  #########
  model = "pop_year;0;5.Cville"
  
  file.cvile <- paste(base, model, 
                      paste(model, ".glmRNP.Rdata", sep = ""), 
                      sep = "/" )
  
  print(file.cvile)
  
  out.glm.cvile.year <- get(load(file.cvile))
  
  
  
  # generate a master index for window analysis
  ### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  
  setkey(out.glm.cvile.year, "chr")
  
  
  
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- out.glm.cvile.year %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  
  wins[,i:=1:dim(wins)[1]]
  
  dim(wins)
  #####  
  glm.out = out.glm.cvile.year
    
  setkey(glm.out, chr, pos)
  head(glm.out)
  
 # wins %<>% filter(chr == "2L")
  
  ### start the summarization process
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind"
  )%do%{
    
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, 
                                    pos=wins[win.i]$start:wins[win.i]$end, 
                                    key="chr,pos")), nomatch=0]
    
    
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]
    
    
    seqSetFilter(genofile, 
                 variant.id=unique(win.tmp$variant.id),
                 sample.id=samps[locality=="VA_ch"][year>= 2016 ]$sampleId)
    
    #obtain AFs 
    af <- seqGetData(genofile, "annotation/format/FREQ")
    f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), 
                        variant.id=seqGetData(genofile, "variant.id"))
    
    #merge AFs with object
    win.tmp <- merge(win.tmp, f.hat, by="variant.id")
    win.tmp[,het:=2*fhat*(1-fhat)]
    
    #thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))[,1]
    
    pr.i <- c(#0.05,
      0.05
    )
    
    #tmpo <- foreach(
    # pr.i=thrs, 
    # .errorhandling="remove", 
    # .combine="rbind")%do%{
    win.tmp %>% 
      filter(!is.na(rnp), pr.i == pr.i ) %>%
      group_by(perm, chr , variable, mod) %>%
      summarise(pos_mean = mean(pos),
                pos_mean = mean(pos),
                pos_min = min(pos),
                pos_max = max(pos),
                win=win.i,
                pr=pr.i,
                rnp.pr=c(mean(rnp<=pr.i)),
                rnp.binom.p=c(binom.test(sum(rnp<=pr.i), 
                                         length(rnp), pr.i)$p.value),
                wZa=sum(het*Z)/(sqrt(sum(het^2))),
                wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                min.p.lrt=min(p_lrt),
                min.rnp=min(rnp),
                nSNPs = n(),
                sum.rnp=sum(rnp<=pr.i),
      ) %>%
      mutate(
        model.pop = model,
        perm_type=ifelse(perm==0, "real","permuted"),
        invName=case_when(
          chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
          chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
          chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
          chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
          TRUE ~ "noInv"
        )) -> win.out
    #}
    #tmpo
  }
  
save(win.out, file = "year.cville.windowAnalysis.Rdata")
load("./year.cville.windowAnalysis.Rdata")

win.out %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.binom.p, 0.00))  %>%
  mutate(metric = "rnvp") -> rnvp

win.out %>%
  mutate(perm_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm",
  )) %>%
  group_by(chr, perm_type, pos_mean) %>%
  summarize(uci = quantile(rnp.wZa.p, 0.00))  %>%
  mutate(metric = "wza") -> wza

rbind(rnvp, wza ) %>%
  ggplot(
    aes(
      x=pos_mean/1e6,
      y=-log10(uci),
      color=perm_type
    )
  ) +
  geom_line(alpha = 0.7) +
  ylab("P-values") +
  xlab("Genome Position (Mb)") +
  theme_bw() +
  scale_color_manual(values = c("grey","red") ) +
  facet_grid(chr~metric, scales = "free") ->
  year.rnp.plot

ggsave(year.rnp.plot, file = "year.rnp.plot.pdf")

#####  

win.out %<>%
  mutate(chr_inv = paste(chr, invName, sep = "_"))

win.out %>%
  group_by(chr,  chr_inv, perm ) %>%
  summarize(me.val = median(rnp.binom.p)) ->
  dat.for.plot

ggplot() +
  geom_violin(
    data=filter(dat.for.plot, perm != 0 ),
    aes(
    x=chr_inv,
    y=-log10(me.val)
  )) +
  geom_point(
    data=filter(dat.for.plot, perm == 0 ),
    aes(
      x=chr_inv,
      y=-log10(me.val)
    ), size = 4, shape = 23, fill = "red"
  )  +
  theme_bw() +
  coord_flip() ->
  violin.perm.dot.real

ggsave(violin.perm.dot.real, file = "violin.perm.dot.real.pdf", w= 4, h =3)



  
  
  
  
  #########
  #########
  #########
  model = "temp.max;2;5.Cville"
  
  file.cvile <- paste(base, model, 
                      paste(model, ".glmRNP.Rdata", sep = ""), 
                      sep = "/" )
  
  print(file.cvile)
  
  out.glm.cvile <- get(load(file.cvile))

### Plot the crossCluster Rnrichment
### data
  file <- "/project/berglandlab/alan/environmental_ombibus_global/crossCluster_enrichment.Rdata"
  load(file)

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

  ggsave(p4, file="./p4_alt_3.pdf", h=7, w=14)

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

p4 <-ggplot(data=m.ag.ag.3[thr==.05][order(-cross.rr.mean)]) +
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
