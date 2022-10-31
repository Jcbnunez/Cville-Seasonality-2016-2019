### libraries
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(viridis)
  library(patchwork)

### data
### scp aob2x@rivanna.hpc.virginia.edu:~/density_analysis_output.Rdata ~/.

load("~/density_analysis_output.Rdata")


rnp.plot <-
ggplot() +
geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
geom_line(data=dens.ag[locality=="VA_ch"][perm!=0][stat=="rnp"],
          aes(x=x, y=y,
              group=interaction(perm, inv),
              linetype=inv),
          color="grey") +
geom_line(data=dens.ag[locality=="VA_ch"][perm==0][stat=="rnp"],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="black", size=1) +
facet_grid(chr~mod) +
theme_bw() +
theme(legend.position="bottom")



chisq.plot <-
ggplot() +
geom_vline(xintercept=qlogis(0.01), linetype="dashed", color="black") +
geom_line(data=dens.ag[locality=="VA_ch"][perm!=0][stat=="chisq"],
          aes(x=x, y=y,
              group=interaction(perm, inv),
              linetype=inv),
          color="grey") +
geom_line(data=dens.ag[locality=="VA_ch"][perm==0][stat=="chisq"],
            aes(x=x, y=y,
                group=interaction(perm, inv),
                linetype=inv),
            color="black", size=1) +
facet_grid(chr~mod) +
theme_bw() +
theme(legend.position="bottom")


dens.ag[locality=="VA_ch"][perm==0][stat=="chisq"][chr=="2L"][mod=="aveTemp+year_factor"]


  mega <-
  rnp.plot + chisq.plot +
  plot_annotation(tag_levels="A")


  ggsave(mega, file="~/densPlot_nested.SNP.pdf", h=10, w=10)
