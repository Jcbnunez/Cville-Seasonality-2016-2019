# ijob -A berglandlab -c30 -p largemem --mem=100G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  #registerDoMC(4)

  registerDoMC(30)

### combine with snp identity
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

  use <- apply(snp.dt[,c("VA_ch"), with=F],
                1, any)
  snp.dt <- snp.dt[use]
  setkey(snp.dt, chr, pos)

### file list
  fl <- list.files("/project/berglandlab/alan/environmental_ombibus_global_permV2/bestAIC/", full.names=T)
  fl <- fl[grepl("v2", fl)]

  best.aic <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[1]
    message(which(fl.i==fl))
    load(fl.i)
    setkey(glm.out.ag, chr, pos)
    tmp <- merge(glm.out.ag, snp.dt)
    return(tmp)
  }

  o.ag <- rbindlist(best.aic)
  setkey(o.ag, chr, pos)

### save
  save(o.ag, file="/project/berglandlab/alan/environmental_ombibus_global_permV2/bestAIC.Rdata")
  load(file="/project/berglandlab/alan/environmental_ombibus_global_permV2/bestAIC.Rdata")

### count numbers
  o.ag.ag <- o.ag[,
                  list(.N),
                  list(perm, cluster, mod, var, chr, inv=invName!="none", stage)]

### summarize
  ### this is how we did it for the first submission
  # o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var, chr, inv),
  #                list(perm, cluster, stage)]

  ### this is how we do it now
    o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var),
                   list(perm, cluster, chr, inv, stage)]

### aggregate
  o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                    prop.perm.mu=median(prop[perm!=0]),
                    prop.perm.lci=quantile(prop[perm!=0], .01),
                    prop.perm.uci=quantile(prop[perm!=0], .99),
                    prop.perm.med=median(prop[perm!=0]),
                    prop.rr=median(log2(prop[perm==0]/prop[perm!=0])),
                    prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
                list(chr, inv=inv, mod, var, cluster, stage)]
  o2.ag[,rr:=prop.real/prop.perm.mu]
  o2.ag[,en:=(prop.real-prop.perm.mu)/prop.perm.mu]


  o2.ag[order(-prop.real)][]
  o2.ag[prop.real>(prop.perm.uci)][order(rr)]
  o2.ag[prop.rr-2*prop.sd>0]
  o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]

  o2.ag[cluster=="1.Europe_W"][order(-rr)][1:10]
  o2.ag[cluster=="2.North_America_W"][order(-rr)][1:10]
  o2.ag[cluster=="2.North_America_E"][order(-rr)][1:10]
  o2.ag[cluster=="3.Europe_E"][order(-rr)][1:10]
  o2.ag[cluster=="5.Cville"][order(-rr)][1:10]
  o2.ag[cluster=="5.Cville"][stage=="stage1"][order(-rr)][1:10]
  o2.ag[cluster=="5.Cville"][stage=="stage2"][order(-rr)][1:10]

  o2.ag[cluster=="5.Cville"][order(-prop.real)][1:10]

### a bit more manipulation
  #o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]
  o2.ag[,sig:=(prop.rr-2*prop.sd)>1 | (prop.rr-2*prop.sd)<1]

  o2.rank <- o2.ag[chr=="2L"][inv==T][cluster=="5.Cville"]
  o2.rank[,modRank:=rank(rr)]
  setkey(o2.rank, mod, var)
  setkey(o2.ag, mod, var)
  o2.ag <- merge(o2.ag, o2.rank[,c("mod", "var", "modRank")])
  o2.ag[,inv:=ifelse(inv, "Inside Inversion", "Outside Inversion")]
  o2.ag[,inv:=factor(inv, levels=c("Inside Inversion", "Outside Inversion"))]

### add in qualifier and also tack in old results.
  o2.ag[,perm_strategy:="new"]
  aic <- o2.ag

### old results
  load(file="/project/berglandlab/alan/environmental_ombibus_global/o2.globalOmnibus.Rdata")
  o2.ag[,perm_strategy:="old"]

### merge
  o2.ag <- rbind(aic, o2.ag, fill=T)

### save
  save(o2.ag, file="/project/berglandlab/alan/environmental_ombibus_global_permV2/bestAIC/bestAIC.v2.Rdata")
