
  library(data.table)
  library(ggplot2)
  library(patchwork)

  inv.dt <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/InversionsMap_hglft_v6_inv_startStop.txt")
  
  setnames(inv.dt, "chrom", "chr")


  load("/project/berglandlab/alan/all_gwas_glm_win3.Rdata")
  
  #gwas.win.o[,grm:=!grepl("nogrms", gwas.pheno)]
  #gwas.win.o.ag <- gwas.win.o[,list(n=sum(pa<.05 & or>1)), list(chr, start, end, grm, glm.perm.n, inv, win.i)]
  gwas.win.o.ag.ag <- gwas.win.o.ag[,list(n=mean(n), sd=sd(n)), list(chr, start, end, grm, perm=I(0!=glm.perm.n), inv)]

  
  ggplot() +
  geom_rect(data=inv.dt[which(inv.dt$invName == "2Lt"),],
            aes(xmin=start/1e6, xmax=stop/1e6, ymin=-1, ymax=20), color="grey", alpha=.5) +
  geom_ribbon(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
                aes(x=I(start/2+end/2)/1e6, 
                    ymax=n+2*sd,
                    ymin=0), 
                color="black", alpha = 0.6) +
  geom_point(data=gwas.win.o.ag.ag[perm==F][grm==F][chr=="2L"],
        aes(x=I(start/2+end/2)/1e6, y=n), color="red") +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
        aes(x=I(start/2+end/2)/1e6, y=n), color="black") +
  #geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
  #      aes(x=I(start/2+end/2)/1e6, y=n-2*sd), color="black", linetype="dashed") +
  facet_grid(grm~chr) +
  theme_bw() -> pheno_plot

ggsave(pheno_plot, file = "pheno_plot.pdf")

  unique(gwas.win.o[glm.perm.n==0][chr=="2L"][pa<.0005][inv==T]$gwas.pheno)
  unique(gwas.win.o[glm.perm.n==0][chr=="2L"][inv==T]$gwas.pheno)
