# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(tidyr)

### import
  fl <- list.files("/sfs/qumulo/qproject/berglandlab/Adam/double_enrichment_slidingwindow", full.names=T)

  gwas.o <- foreach(fl.i=fl)%dopar%{

    message(fl.i)
    load(fl.i)
    gwas.win.o[,perm:=tstrsplit(glm.perm, "\\.")[[3]]]
    gwas.win.o[,perm:=tstrsplit(perm, "_")[[3]]%>%as.numeric]

    pa <- gwas.win.o[,list(fet.pa=p.adjust(fet.p), prop.pa=p.adjust(prop.p), win.i), list(gwas.pheno, glm.perm)]
    setkey(gwas.win.o, gwas.pheno, glm.perm, win.i)
    setkey(pa, gwas.pheno, glm.perm, win.i)
    gwas.win.o <- merge(gwas.win.o, pa)


    #gwas.win.o
    return(gwas.win.o)
  }

  gwas.o <- rbindlist(gwas.o)
  gwas.o[,grm:=!grepl("nogrms", gwas.pheno)]

  gwas.o.ag <- gwas.o[,list(nSig.en=sum(fet.pa[or>1]<.05), nSig.dep=sum(fet.pa[or<1]<.05), nSig.prop=sum(prop.pa<.05)), list(grm, chr, inv, gwas.pheno, perm)]

  save(gwas.o.ag, file="/project/berglandlab/alan/pheno.glm.chr.Rdata")

  summary(lm(nSig.dep~as.factor(perm), gwas.o.ag[grm==T]))

  setkey(gwas.o, grm)
  tab <- table(gwas.o[J(T)]$fet.pa<.05, gwas.o[J(T)]$perm, gwas.o[J(T)]$chr)

  save(gwas.o, file="/project/berglandlab/alan/pheno.glm.Rdata")


  fet.fun <- function(mat){
    fet <- fisher.test(mat)
    fet$p.value
  }

fet.p=fet.fun(mat=matrix(c(TT, TF, FT, FF), nrow=2, byrow=T)),

  gwas.o.ag <- gwas.o[win.i%%2==0,list(TT=sum(TT), TF=sum(TF), FT=sum(FT), FF=sum(FF), st.T=sum(st.T), st.F=sum(st.F)), list(chr, inv, gwas.pheno, perm, grms=!grepl("nogrms", gwas.pheno))]
  gwas.o.ag.p <- gwas.o.ag[,list(
                                 fet.or=(TT/TF)/(FT/FF),
                                 st.pr=(st.T)/(st.T+st.F)),
                            list(chr, inv, gwas.pheno, perm, grms)]

  gwas.o.ag.ag <- gwas.o.ag.p[,list(fet.pr=mean(fet.or[perm==0]>fet.or[perm!=0], na.rm=T), fet.or=fet.or[perm==0],
                                    st.pr=mean(st.pr[perm==0]>st.pr[perm!=0], na.rm=T), st.prop=st.pr[perm==0]),
                               list(chr, inv, gwas.pheno, grms)]


  save(gwas.o.ag.ag, file="~/gwas_glm_summary.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/gwas_glm_summary.Rdata ~/.


  load("~/gwas_glm_summary.Rdata")

  ggplot(data=gwas.o.ag.ag[fet.or<5][grms==T]) +
  geom_point(data=gwas.o.ag.ag[fet.or<5][grms==T], aes(x=log2(fet.or), y=st.prop))+
  geom_point(data=gwas.o.ag.ag[fet.or<5][fet.pr>.95][grms==T], aes(x=log2(fet.or), y=st.prop), color="red", size=.5) +
  facet_grid(inv~chr)



scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/pheno.glm.chr.Rdata ~/.

gwas.o.ag.ag <- gwas.o.ag[,list(nSig.en.mu=mean(nSig.en[perm=!0], na.rm=T),
                                nSig.en.sd=sd(nSig.en[perm=!0], na.rm=T),
                            nSig.en=nSig.en[perm=!0]),
                          list(chr, inv, gwas.pheno, grm)]

ggplot(data=gwas.o.ag.ag) +
geom_point(aes(x=as.numeric(as.factor(gwas.pheno)), y=nSig.en.mu)) +
geom_linerange(aes(x=as.numeric(as.factor(gwas.pheno)), ymin=nSig.en.mu-2*nSig.en.sd, ymax=nSig.en.mu+2*nSig.en.sd)) +
geom_point(aes(x=as.numeric(as.factor(gwas.pheno)), y=nSig.en), color="red") +

facet_grid(inv~chr)


table(gwas.o.ag.ag[grms==T]$fet.pr>.99, gwas.o.ag.ag[grms==T]$chr, gwas.o.ag.ag[grms==T]$inv)
table(gwas.o.ag.ag[grms==F]$fet.pr>.99, gwas.o.ag.ag[grms==F]$chr, gwas.o.ag.ag[grms==F]$inv)

### copy to load and load
### scp aob2x@rivanna.hpc.virginia.edu:~/all_gwas_glm2.Rdata ~/.

  library(data.table)
  library(ggplot2)
  library(patchwork)

  load("~/all_gwas_glm2.Rdata")
  gwas.o[,grm:=!grepl("nogrms", gwas.pheno)]

### basic plot
  p1 <-
  ggplot(data=gwas.o[glm.perm==0][grm==T], aes(x=log10(thr), y=log2(or), group=interaction(inv, gwas.pheno), color=inv)) +
  geom_line() +
  facet_grid(chr~glm.mod)

  p2 <-
  ggplot(data=gwas.o[glm.perm==0][grm==F], aes(x=log10(thr), y=log2(or), group=interaction(inv, gwas.pheno), color=inv)) +
  geom_line() +
  facet_grid(chr~glm.mod)

  mega <- p1 + p2

  ggsave(mega, file="~/gwas.glm.pdf", h=12, w=20)

### simpler plot @ fixed threshold
  or.plot <-
  ggplot() +
  geom_density(data=gwas.o[glm.perm!=0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>50],
              aes(x=log2(or), group=glm.perm, color=as.factor(glm.perm==0))) +
  geom_density(data=gwas.o[glm.perm==0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>50],
              aes(x=log2(or), group=glm.perm, color=as.factor(glm.perm==0))) +
  facet_grid(chr~inv)

  pr.plot <-
  ggplot() +
  geom_density(data=gwas.o[glm.perm!=0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>50],
              aes(x=prop, group=glm.perm, color=as.factor(glm.perm==0))) +
  geom_density(data=gwas.o[glm.perm==0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>50],
              aes(x=prop, group=glm.perm, color=as.factor(glm.perm==0))) +
  facet_grid(chr~inv)

  mega <-
  or.plot + pr.plot
  ggsave(mega, file="~/gwas.glm.pdf", h=12, w=20)


### OR prop biplot
  gwas.o.ag <- gwas.o[thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>20][,
                      list(or.mean=mean(log2(or)), prop.mean=mean(prop),
                           or.sd=sd(log2(or)), prop.sd=sd(prop)),
                      list(perm=as.factor(glm.perm>0), gwas.pheno, chr, inv)]

  ggplot() +
  geom_point(data=gwas.o[glm.perm!=0][thr==0.01][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>10],
              aes(x=log2(or), y=prop), color="black", alpha=.05) +
  geom_point(data=gwas.o[glm.perm==0][thr==0.01][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>10],
              aes(x=log2(or), y=prop), color="red") +
  facet_grid(chr~inv)


### empirical probabilities
  gwas.o.ag <- gwas.o[st.T+st.F>10,list(pr.or=mean(or[glm.perm==0]>or[glm.perm!=0]),
                            pr.prop=mean(prop[glm.perm==0]>prop[glm.perm!=0])),
                       list(gwas.pheno, glm.mod, thr, chr, inv, grm)]


 ggplot() +
 geom_point(data=gwas.o.ag[thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][],
             aes(x=pr.or, y=pr.prop)) +
 facet_grid(chr~inv) +
 geom_vline(xintercept=.9) +
 geom_hline(yintercept=c(.1, .9))


gwas.o.ag[thr==0.01][grm==F][glm.mod=="aveTemp+year_factor"][pr.or>.9 & (pr.prop>.9)][chr=="2L"][inv==T][order(pr.prop)][,c("gwas.pheno", "pr.or", "pr.prop"), with=F]
gwas.o.ag[thr==0.01][grm==F][glm.mod=="aveTemp+year_factor"][pr.or>.9 & (pr.prop<.1)][chr=="2L"][inv==T][order(pr.prop)][,c("gwas.pheno", "pr.or", "pr.prop"), with=F]



gwas.o[gwas.pheno=="StarvationResistance_standard_male.nogrms.original"][glm.mod=="aveTemp+year_factor"][thr==0.05][chr=="2L"][inv==T][,
list(pr=mean(prop[glm.perm==0]>prop[glm.perm!=0]))]



gwas.o[glm.perm==0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>20][log2(or)>.5]



table(gwas.o[grm==F][chr=="2L"][thr==0.05][glm.mod=="aveTemp+year_factor"]$fet.p<0.00001,
      gwas.o[grm==F][chr=="2L"][thr==0.05][glm.mod=="aveTemp+year_factor"]$glm.perm,
      gwas.o[grm==F][chr=="2L"][thr==0.05][glm.mod=="aveTemp+year_factor"]$inv)


  ggplot() +
  geom_point(data=gwas.o.ag,
              aes(x=or.mean, y=prop.mean, color=perm)) +
  facet_grid(chr~inv)
