# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(tidyr)
  library(doMC)
  registerDoMC(5)

### import
  fl <- list.files("/scratch/aob2x/gwas_glm_merge_window/", full.names=T)

  gwas.win.o <- foreach(fl.i=fl)%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    gwas.win.o[,glm.perm.n:=tstrsplit(glm.perm, "/")%>%last]
    gwas.win.o[,glm.perm.n:=tstrsplit(glm.perm.n, "\\.")[[3]]]
    gwas.win.o[,glm.perm.n:=tstrsplit(glm.perm.n, "_")[[3]] %>% as.numeric]
    gwas.win.o[,pa:=p.adjust(gwas.win.o$fet.p)]

    return(gwas.win.o)
  }
  gwas.win.o <- rbindlist(gwas.win.o)


  save(gwas.win.o, file="~/all_gwas_glm_win.Rdata")

### noodling

  table(gwas.win.o$pa<.05, gwas.win.o$glm.perm.n, gwas.win.o$chr)

  table(gwas.win.o[chr=="2L"]$pa<.05, gwas.win.o[chr=="2L"]$glm.perm.n, gwas.win.o[chr=="2L"]$inv)

  table(gwas.win.o[chr=="2L"][glm.perm.n==0][pa<.05][inv==T]$gwas.pheno)

### copy to load and load
### scp aob2x@rivanna.hpc.virginia.edu:~/all_gwas_glm_win.Rdata ~/.

  library(data.table)
  library(ggplot2)
  library(patchwork)


  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")


  load("~/all_gwas_glm_win.Rdata")
  gwas.win.o[,grm:=!grepl("nogrms", gwas.pheno)]
  #gwas.win.o[,pa:=p.adjust(gwas.win.o$fet.p)]

  gwas.win.o.ag <- gwas.win.o[,list(n=sum(pa<.05 & or>1)), list(chr, start, end, grm, glm.perm.n, inv, win.i)]
  gwas.win.o.ag.ag <- gwas.win.o.ag[,list(n=median(n), sd=sd(n)), list(chr, start, end, grm, perm=I(0!=glm.perm.n), inv)]

  ggplot() +
  geom_rect(data=inv.dt, aes(xmin=start/1e6, xmax=stop/1e6, ymin=-1, ymax=20), color="grey", alpha=.5) +
  geom_point(data=gwas.win.o.ag.ag[perm==F][grm==F],
        aes(x=I(start/2+end/2)/1e6, y=n), color="red") +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F],
        aes(x=I(start/2+end/2)/1e6, y=n), color="black") +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F],
        aes(x=I(start/2+end/2)/1e6, y=n-2*sd), color="black", linetype="dashed") +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F],
        aes(x=I(start/2+end/2)/1e6, y=n+2*sd), color="black", linetype="dashed") +
  facet_grid(grm~chr) +
  theme_bw()



  +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName))


  ggplot(data=gwas.win.o[glm.perm.n==0][pa<.05],
        aes(x=I(start/2+end/2)/1e6, y=-log10(pa))) +
  geom_point() +
  facet_grid(~chr)+
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName))
































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
  ggplot() +
  geom_point(data=gwas.o[glm.perm!=0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>20],
              aes(x=log2(or), y=prop), color="black", alpha=.5) +
  geom_point(data=gwas.o[glm.perm==0][thr==0.05][grm==F][glm.mod=="aveTemp+year_factor"][st.T+st.F>20],
              aes(x=log2(or), y=prop), color="red") +

  facet_grid(chr~inv)
