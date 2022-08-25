scp aob2x@rivanna.hpc.virginia.edu:~/perm.dest_core20.Rdata ~/.

library(data.table)
library(ggplot2)
library(patchwork)

load("~/perm.dest_core20.Rdata")
o.ag <- o.dest_core20[,list(or.mu=mean(log2(or)), or.lci=quantile(log2(or), .025), or.uci=quantile(log2(or), .975),
                            pr.mu=mean(st.pr), pr.lci=quantile(st.pr, .025), pr.uci=quantile(st.pr, .975)),
                       list(perm=I(perm!=0), chr, inv, thr, mod)]



p1 <- ggplot(data=o.ag[thr==0.05][mod=="core20_DEST"]) +
geom_linerange(data=o.ag[thr==0.05][mod=="core20_DEST"][perm==T],
          aes(x=interaction(chr, inv), ymin=or.lci, ymax=or.uci, color=as.factor(perm))) +
geom_point(aes(x=interaction(chr, inv), y=or.mu, color=as.factor(perm)))



p2 <- ggplot(data=o.ag[thr==0.05][mod=="core20_DEST"]) +
geom_linerange(data=o.ag[thr==0.05][mod=="core20_DEST"][perm==T],
          aes(x=interaction(chr, inv), ymin=pr.lci, ymax=pr.uci, color=as.factor(perm))) +
geom_point(aes(x=interaction(chr, inv), y=pr.mu, color=as.factor(perm)))


p1+p2
