
  library(data.table)
  library(ggplot2)
  library(patchwork)


  inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")

  #### raw
    load("~/all_gwas_glm_win2.Rdata")
    gwas.win.o[,grm:=!grepl("nogrms", gwas.pheno)]
    gwas.win.o[,nSNPs:=TT+TF+FT+FF]

    gwas.win.o.ag <- gwas.win.o[nSNPs>250, list(n=sum(pa<.05), terms=paste(gwas.pheno[pa<.05], collapse=";")),
                                list(chr, start, end, grm, glm.perm.n, inv, win.i)]

    gwas.win.o.ag.ag <- gwas.win.o.ag[,list(n=mean(n), sd=sd(n), terms=paste(terms, collapse="|")),
                                            list(chr, start, end, grm, win.i, perm=I(0!=glm.perm.n), inv)]


    unlist(tstrsplit(gwas.win.o.ag.ag[n>25]$terms, ";"))

  ### summarized
  load("~/all_gwas_glm_win3.Rdata")


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


####


gwas.win.o.ag.ag[perm==F][n>15]






### by trait
  table(
    gwas.win.o[nSNPs>250 & grm==F][glm.perm.n==0]$pa<.05,
    gwas.win.o[nSNPs>250 & grm==F][glm.perm.n==0]$or<1)


  gwas.win.o.ag2 <- gwas.win.o[nSNPs>250 & grm==F,
                              list(n=sum(pa<.05 & or>1)),
                              list(chr, gwas.pheno, grm, win.i, glm.perm.n, inv)]


  tmp <- gwas.win.o.ag2[glm.perm.n!=0][chr=="2L"][inv==F][gwas.pheno=="ChillComaRecoveryTime_standard_female.nogrms.original"]
  e <- gwas.win.o.ag2[,list(n.obs=mean(n[glm.perm.n==0]>=1), n.exp=mean(n[glm.perm.n!=0]>1)),
                                      list(chr, gwas.pheno, grm, inv)]


  gwas.win.o.ag2.ag.ag <- gwas.win.o.ag2.ag[,list(.N, prop=mean(n>0)), list(perm, chr, inv, grm)]


ggplot(gwas.win.o.ag2.ag.ag) +
geom_point(aes(x=perm, y=prop, color=as.factor(perm==0))) +
facet_grid(chr~inv) +




  table(gwas.win.o.ag2.ag$n>0, gwas.win.o.ag2.ag$perm, gwas.win.o.ag2.ag$chr, gwas.win.o.ag2.ag$inv)


  p1 <- ggplot(data=gwas.win.o.ag2.ag[n>.01]) +
  geom_point(aes(x=n, y=gwas.pheno, color=as.factor(perm))) +
  facet_grid(inv~chr) +
  theme(legend.position = "none", axis.text.y=element_text(size=4))

  ggsave(p1, file="~/pheno_enrichr.pdf")


  geom_linerange(aes(xmin=n-2*sd, xmax=n+2*sd, y=gwas.pheno, color=perm)) +

  tmp <- gwas.win.o.ag2.ag[perm==F & chr=="2L" & inv==T,list(n=n), list(gwas.pheno)]
  tmp[,r:=rank(n, ties="random")]
  gwas.win.o.ag2.ag[,gwas.pheno:=factor(gwas.pheno, levels=tmp[order(r)]$gwas.pheno)]

  ggplot() +
  geom_errorbar(data=gwas.win.o.ag2.ag[inv==T][chr=="2L"][perm==T],
            aes(xmin=n-2*sd, xmax=n+2*sd, y=gwas.pheno, group=as.factor(perm), color=as.factor(perm))) +
  geom_point(data=gwas.win.o.ag2.ag[inv==T][chr=="2L"],
            aes(y=gwas.pheno, x=n, group=as.factor(perm), color=as.factor(perm)), size=.5) +
  theme(legend.position = "none", axis.text.y=element_text(size=1))



  gwas.win.o.ag2.ag.ag <- gwas.win.o.ag2.ag[,list(delta=n[perm==F] - n[perm==T]),
                                      list(chr, gwas.pheno, grm, inv)]


  ggplot(data=gwas.win.o.ag2.ag.ag, aes(delta)) +
  geom_histogram() +
  facet_grid(inv~chr)

gwas.win.o.ag2.ag[perm==F][n>.01][chr=="2L"][inv==T]


gwas.win.o.ag2.ag.ag[chr=="2L"][inv==T][delta>.02]

setkey(gwas.win.o, win.i)
  focalWindows <- unique(gwas.win.o.ag.ag[perm==F][chr=="2L"][n>12][inv==T]$win.i)
  gwas.win.o[glm.perm.n==0 & grm==F][J(focalWindows)][pa<.05][or>1]$gwas.pheno %>% unique

  focalWindows <- unique(gwas.win.o.ag.ag[perm==F][chr=="2L"][n>9][inv==T][win.i>75]$win.i)
  gwas.win.o[glm.perm.n==0 & grm==F][J(focalWindows)][pa<.05][or>1]$gwas.pheno %>% unique


  focalWindows <- unique(gwas.win.o.ag.ag[perm==F][chr=="2R"][n>4]$win.i)
  gwas.win.o[glm.perm.n==0 & grm==F][J(focalWindows)][pa<.05][or>1]$gwas.pheno %>% unique


  ggplot(data=gwas.win.o[grm==T]) +
  geom_density(aes(-log10(fet.p), group=glm.perm.n, color=I(glm.perm.n==0))) +
  facet_grid(inv~chr)



  unique(gwas.win.o[glm.perm.n==0][chr=="2L"][inv==T]$gwas.pheno)
