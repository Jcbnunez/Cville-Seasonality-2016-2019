# ijob -A berglandlab_standard -c20 -p standard --mem=20G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### combine with snp identity
  load("~/Overwintering_18_19/SNP_filtering/snp_dt_25percMissing.Rdata")
  setkey(snp.dt, id)

  use <- apply(snp.dt[,c("VA_ch"), with=F],
                1, any)
  snp.dt <- snp.dt[use]
  setkey(snp.dt, id)


### split output and get summaries
  fl <- list.files("/project/berglandlab/alan/environmental_ombibus/rawData",full.names=T)
  #fl <- paste("/scratch/aob2x/environmental_ombibus/job", c(1:5000), ".Rdata", sep="")

  mainDir <- "/scratch/aob2x/environmental_ombibus"
  dir.create(file.path(mainDir), showWarnings = FALSE)

  fl.i <- fl[1]
  message(fl.i)
  load(fl.i)
  mod_var <- unique(paste(glm.out$variable, glm.out$mod, sep="_"))


  o.ag <- foreach(fl.i=fl, .errorhandling="remove")%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    glm.out <- merge(glm.out, snp.dt, by.x="variant.id", by.y="id", all.x=T)
    setkey(glm.out, variable, mod)


    glm.out.ag <- glm.out[,
                        list(var=variable[which.min(AIC)], mod=mod[which.min(AIC)],
                             p_lrt=p_lrt[which.min(AIC)],
                             minAIC=min(AIC), yearAIC=AIC[variable=="year"]),
                        list(variant.id, perm)]
    glm.out.ag

    return(glm.out.ag)

  }
  o.ag <- rbindlist(o.ag)

  o.ag <- merge(o.ag, snp.dt, by.x="variant.id", by.y="id")
  o.ag.ag <- o.ag[,list(delta=mean(minAIC - yearAIC)), list(chr, inv=invName!="none", mod, var, perm)]


  save(o.ag, file="/project/berglandlab/alan/environmental_ombibus/bestAIC.Rdata")
  save(mod_var, file="/project/berglandlab/alan/environmental_ombibus/mod_var.Rdata")






































### rank normalization test
  foreach(mvi=mod_var)%do%{
    mvi.dir <-  paste(mainDir, mvi, sep="/")
    fl <- list.files(mvi.dir, full.names=T)
    glm.out <- foreach(fl.i=fl[grepl("job")])%dopar%{
      message(fl.i)
      #fl.i <- fl[1]
      tmp <- fread(fl.i)
      tmp
    }
    glm.out <- rbindlist(glm.out)

    glm.out.ag <- glm.out[,list(rnp=rank(p_lrt)/length(p_lrt), chr, pos, inv=invName!="none"), list(perm)]

    thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-2)))[,1]
    #thrs <- c(.05, .01)

    o.rnp.ag <- foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%do%{
      message(thr.i)
      #thr.i <- thrs[1]
      o.rnp.ag <- glm.out.ag[,list(pr=mean(rnp<thr.i), .N,
                                   thr=thr.i),
                        list(chr, inv, perm)]
      o.rnp.ag
    }
    o.rnp.ag[,mod:=glm.out[1]$mod]
    o.rnp.ag[,variable:=glm.out[1]$variable]
    save(o.rnp.ag, file=paste(mvi.dir, "rnp_summary.Rdata", sep="/"))

  }
  fl <- list.files("/scratch/aob2x/dest_glm_final_nested_qb_alternatives_AIC",full.names=T)























### summarize - how frequently are is each model the Best model
  o.ag.ag <- o.ag[,
                  list(N=length(variant.id)),
                  list(perm, mod, var, chr, inv=invName!="none")]

  o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var),
                      list(perm , chr, inv)]

  o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                    prop.perm.mu=mean(prop[perm!=0]),
                    prop.perm.lci=quantile(prop[perm!=0], .01),
                    prop.perm.uci=quantile(prop[perm!=0], .99),
                    prop.perm.med=median(prop[perm!=0]),
                    prop.rr=mean(log2(prop[perm==0]/prop[perm!=0])),
                    prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
                list(chr, inv, mod, var)]
  o2.ag[,rr:=prop.real/prop.perm.mu]

  o2.ag[order(-prop.real)][]
  o2.ag[prop.real>(prop.perm.uci)][order(rr)]
  o2.ag[prop.rr-2*prop.sd>0]
  o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]

### save
  save(o2.ag, file="~/environmental_omnibus_summary.Rdata")































### what is the average deltaAIC for the best model?

  o3.ag <- o.ag[,list(mean_deltaAIC_year=mean(deltaAIC_year),
                         sd_deltaAIC_year=sd(deltaAIC_year)),
                    list(chr, inv=I(invName!="none"), mod, perm)]

  o3.ag[order(mean_deltaAIC_year)]
  #o3.ag[perm==F][mod%in%c(3,4)]
  #o3.ag[mod==4][chr=="2L"]

### RNP for each mod
  o.rnp <- o.ag[N==0 & cm_mb>0 & !is.na(cm_mb),
                  list(p.lrt.modF4, p.lrt.modD3,
                      rank_modF4=rank(p.lrt.modF4)/length(p.lrt.modF4),
                      rank_modD3=rank(p.lrt.modD3)/length(p.lrt.modD3),
                      variant.id=variant.id, chr, inv=I(invName!="none"), mod), list(perm)]

  thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-2)))[,1]
  #thrs <- c(.05, .01)
  o.rnp.ag <- foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%do%{
    message(thr.i)
    #thr.i <- thrs[1]
    o.rnp.ag <- o.rnp[,list(pr_modF4=mean(rank_modF4<thr.i),
                            pr_modD3=mean(rank_modD3<thr.i), thr=thr.i),
                      list(chr, inv, perm)]
    o.rnp.ag
  }


  o.rnp.ag.ag <- o.rnp.ag[,list(rnp.pr=c(pr_modF4[perm==0], pr_modD3[perm==0]),
                                rnp.pr.perm=c(mean(pr_modF4[perm!=0]), mean(pr_modD3[perm!=0])),
                                mod=LETTERS[c(5,4)]),
                            list(chr, inv, thr)]



    summary(lm(log10(p.lrt.modF4) ~ log10(p.lrt.modD3), o.rnp[perm==0]))

    table(o.rnp[perm==0]$rank_modF4<.01, o.rnp[perm==0]$rank_modD5<.01)

### save
  save(o2.ag, o3.ag, o.rnp.ag, o.rnp.ag.ag, file="~/o2_modSelection.Rdata")

### save
  glm.out.ms <- o.ag[perm==0]
  save(glm.out.ms, file="/scratch/aob2x/modSel_out.Rdata")


### how much more significant is modF4 than modE3?
  prop.table(table(
  o.ag[perm==0]$p.lrt.modF4<.0005,
  o.ag[perm==0]$chr,
  o.ag[perm==0]$inv!="none"))

  prop.table(table(
  o.ag[perm==0]$p.lrt.modD3<.0005,
  o.ag[perm==0]$chr,
  o.ag[perm==0]$inv!="none"))

  #save(o.ag, file="~/modelSelection.Rdata")


  table(
  sign(o.ag[perm==0][p.lrt.modF4<.05]$b_temp.modF4) ==
  sign(o.ag[perm==0][p.lrt.modF4<.05]$b_temp.modD3))







  ### load core20 GLM object
    load(file="/project/berglandlab/alan/core20glm.Rdata")

  ### merge
    m1 <- merge(glm.out.ms, core20.glm[mod=="season+locality_factor"][set=="machado"],  by="variant.id")
    m2 <- merge(glm.out.ms, core20.glm[mod=="season+locality_factor"][set=="no_va"],    by="variant.id")
    m1.beta <- m1[,list(seas.beta=tstrsplit(b, ";")%>%last%>%as.numeric), list(variant.id)]
    m2.beta <- m2[,list(seas.beta=tstrsplit(b, ";")%>%last%>%as.numeric), list(variant.id)]

    m1 <- merge(m1, m1.beta, by="variant.id")
    m2 <- merge(m2, m2.beta, by="variant.id")

    m1[!is.na(p.lrt.modF4) & !is.na(p.lrt) ,thermal.rank := rank(p.lrt.modF4)/length(p.lrt.modF4)]
    m1[!is.na(p.lrt.modF4) & !is.na(p.lrt) ,season.rank := rank(p.lrt)/length(p.lrt)]


    m2[!is.na(p.lrt.modF4) & !is.na(p.lrt) ,thermal.rank := rank(p.lrt.modF4)/length(p.lrt.modF4)]
    m2[!is.na(p.lrt.modF4) & !is.na(p.lrt) ,season.rank := rank(p.lrt)/length(p.lrt)]

    thr <- .01
    table(m1$thermal.rank<thr, m1$season.rank<thr)%>%fisher.test()
    table(m2$thermal.rank<thr, m2$season.rank<thr)%>%fisher.test()

    m1[,inv:=invName!="none"]
    m2[,inv:=invName!="none"]

    thrs <- c(0.01, 0.05)

    o <- foreach(chr.i=unique(m1$chr), .errorhandling="pass")%dopar%{
      foreach(inv.i=c(T,F), .errorhandling="pass")%do%{
        foreach(thr.i=thrs, .errorhandling="pass")%do%{
          # chr.i <- "2L"; inv.i <- T; thr.i<-0.05
            message(paste(chr.i, inv.i, thr.i, sep=" / "))

            ### Machado set, enrichment
              tab <- table(m1[chr==chr.i][inv==inv.i]$thermal.rank  <thr.i,
                           m1[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

              fet <- fisher.test(tab)

           ### Machado set, sign test
             st.T <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$b_temp.modF4) ==
                          sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
             st.F <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$b_temp.modF4) !=
                          sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

             bt <- binom.test(st.T, st.T+st.F, .5)


             tmp1 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i, perm=job,
                        mod="core20_Machado", set="aveTemp_yearFactor",
                        or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                        st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])


           ### DEST set, enrichment
             tab <- table(m2[chr==chr.i][inv==inv.i]$thermal.rank  <thr.i,
                          m2[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

             fet <- fisher.test(tab)

           ### DEST set, sign test
             st.T <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$b_temp.modF4) ==
                          sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
             st.F <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$b_temp.modF4) !=
                          sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

             bt <- binom.test(st.T, st.T+st.F, .5)


             tmp2 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i, perm=job,
                        mod="core20_DEST", set="aveTemp_yearFactor",
                        or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                        st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])

          ### return
            rbind(tmp1, tmp2, fill=T)

        }
      }
    }




    setkey(gtw.temp, chr, pos)
    setkey(glm.out.ms, chr, pos)
    gtw.temp <- merge(gtw.temp, glm.out.ms[, c("chr", "pos", "p.lrt.modF4", "b_temp.modF4"),with=F])
    gtw.temp[,rnp.modF4:=rank(p.lrt.modF4)/length(p.lrt.modF4)]



    table(gtw.temp[chr=="2L"]$rnp.modF4<.005, gtw.temp[chr=="2L"]$rnp.temp.NorthAmerica<.005)%>%fisher.test()

    table(gtw.temp[chr=="2L"][inv==F]$rnp.modF4<.005, gtw.temp[chr=="2L"][inv==F]$rnp.lat.NorthAmerica<.005)%>%fisher.test()

    table(sign(gtw.temp[chr=="2L" & rnp.modF4<.05 & rnp.lat.NorthAmerica<.05]$lat.coef_NorthAmerica) ==
          sign(gtw.temp[chr=="2L" & rnp.modF4<.05 & rnp.lat.NorthAmerica<.05]$b_temp.modF4))%>%prop.test()
