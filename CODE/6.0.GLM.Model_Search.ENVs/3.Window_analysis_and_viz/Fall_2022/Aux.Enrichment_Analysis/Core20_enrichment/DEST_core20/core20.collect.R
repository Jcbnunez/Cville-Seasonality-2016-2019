# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)

### import
  fl <- list.files("/scratch/aob2x/dest_core20", full.names=T)

  core20.glm <- foreach(fl.i=fl)%do%{
    message(fl.i)
    o <- fread(fl.i)
    return(o)
  }

  core20.glm <- rbindlist(core20.glm)

  save(core20.glm, file="/project/berglandlab/alan/core20glm.Rdata")






  job <- 0
  glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
  load(glm.fn)


  m1 <- merge(glm.out[mod=="aveTemp+year_factor"], core20.glm[mod=="season+locality_factor"][set=="machado"], by.x="id", by.y="variant.id")
  m2 <- merge(glm.out[mod=="aveTemp+year_factor"], core20.glm[mod=="season+locality_factor"][set=="no_va"], by.x="id", by.y="variant.id")


  rm(glm.out, core20.glm)

  load("/project/berglandlab/jcbnunez/Shared_w_Alan/haplo_tags_SNPids.Rdata")
  haplo_tags_SNPids_and_inv <- as.data.table(haplo_tags_SNPids_and_inv)
  haplo_tags_SNPids_and_inv[,pos:=as.numeric(pos)]
  setkey(m1, chr, pos)
  setkey(m2, chr, pos)

  setkey(haplo_tags_SNPids_and_inv, chr, pos)

  m1 <- merge(m1, haplo_tags_SNPids_and_inv, all.x=T)
  m2 <- merge(m2, haplo_tags_SNPids_and_inv, all.x=T)

  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]


  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]

  m2[is.na(win), win:="outside"]

  tab1 <- table(m1$p.lrt.y<.005, !is.na(m1$class))
  fisher.test(tab1)

  tab2 <- table(m2$p.lrt.y<.05, m2$win)
  fisher.test(tab2)





  tab2 <- table( m1$thermal.rank<.01,
                 m1$season.rank<.01,
                 m1$chr)
  fisher.test(tab2[,,1])


  tab2 <- table( m1[chr=="2L"]$thermal.rank<.01,
                 m1[chr=="2L"]$season.rank <.01,
                 m1[chr=="2L"]$invName)
  fisher.test(tab2[,,1])
  fisher.test(tab2[,,2])

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
