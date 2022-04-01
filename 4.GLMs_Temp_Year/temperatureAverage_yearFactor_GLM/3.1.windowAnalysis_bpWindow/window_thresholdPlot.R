scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest_glm_morePerms_nested_qb/windowAnalysis_lotsPR/WZA_window.dt.VA_ch_0.Rdata ~/.


library(ggplot2)
library(data.table)
library(patchwork)
load("~/WZA_window.dt.VA_ch_0.Rdata")

yf <- ggplot(data=win.out[mod=="year_factor"][nSNPs>50],
            aes(x=log10(pr), y=log10(rnp.pr), color=as.factor(rbinom.p<.05/2000), group=i), alpha=.5) +
geom_line() +
geom_abline(aes(slope=1, intercept=0), color="red") +
facet_grid(chr.x~as.factor(invName!="none")) + theme(legend.position = "none") +
geom_vline(aes(xintercept=log10(0.05))) +
ggtitle("year factor")


at <- ggplot(data=win.out[mod=="aveTemp+year_factor"][nSNPs>50],
            aes(x=log10(pr), y=log10(rnp.pr), color=as.factor(rbinom.p<.05/2000), group=i), alpha=.5) +
geom_line() +
geom_abline(aes(slope=1, intercept=0), color="red") +
facet_grid(chr.x~as.factor(invName!="none")) + theme(legend.position = "none") +
geom_vline(aes(xintercept=log10(0.05))) +
ggtitle("aveTemp")

mega <- yf + at
ggsave(mega, file="~/mega_pr.pdf")
