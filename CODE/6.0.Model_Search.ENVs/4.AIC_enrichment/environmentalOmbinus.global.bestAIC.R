# ijob -A berglandlab_standard -c20 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### combine with snp identity
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

### split output and get summaries
  fl <- list.files("/scratch/aob2x/environmental_ombibus_global/",full.names=T)

  fl.i="/project/berglandlab/alan/environmental_ombibus_global/temp.propMax;11;2.North_America_W"  
  
  o.ag <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[9640]
    message(fl.i)
    load(fl.i)

    glm.out.ag <- glm.out[,
                        list(var=variable[which.min(AIC)], mod=mod[which.min(AIC)],
                             p_lrt=p_lrt[which.min(AIC)],
                             minAIC=min(AIC), yearAIC=AIC[variable=="pop_year"]),
                        list(variant.id, perm, cluster)]

    ### rename new East coast set
      if(any(glm.out.ag$cluster=="2.North_America_Mid")) {
        glm.out.ag[cluster=="2.North_America_E", cluster:="2.North_America_I95"]
      }

    #glm.out.ag <- merge(glm.out.ag, snp.dt, by.x="variant.id", by.y="id", all.x=T)
    glm.out.ag


  }
  o.ag <- rbindlist(o.ag)

  o.ag <- merge(o.ag, snp.dt, by.x="variant.id", by.y="id")

### save raw output
  save(o.ag, file="/project/berglandlab/alan/environmental_ombibus_global/bestAIC.global.Rdata")

### summarize: how many times is each model the best model?
  o.ag.ag <- o.ag[,
                  list(.N),
                  list(perm, cluster, mod, var, chr, inv=invName!="none")]

  o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var, chr, inv),
                 list(perm, cluster)]

  o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                    prop.perm.mu=median(prop[perm!=0]),
                    prop.perm.lci=quantile(prop[perm!=0], .01),
                    prop.perm.uci=quantile(prop[perm!=0], .99),
                    prop.perm.med=median(prop[perm!=0]),
                    prop.rr=median(log2(prop[perm==0]/prop[perm!=0])),
                    prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
                list(chr, inv=inv, mod, var, cluster)]
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

### a bit more manipulation
  o2.ag[,sig:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]

  o2.rank <- o2.ag[chr=="2L"][inv==T][cluster=="5.Cville"]
  o2.rank[,modRank:=rank(rr)]
  setkey(o2.rank, mod, var)
  setkey(o2.ag, mod, var)
  o2.ag <- merge(o2.ag, o2.rank[,c("mod", "var", "modRank")])
  o2.ag[,inv:=ifelse(inv, "Inside Inversion", "Outside Inversion")]
  o2.ag[,inv:=factor(inv, levels=c("Inside Inversion", "Outside Inversion"))]

  save(o2.ag, file="/project/berglandlab/alan/environmental_ombibus_global/o2.globalOmnibus.Rdata")
