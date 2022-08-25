### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load data
  # scp aob2x@rivanna.hpc.virginia.edu:~/dest_core20.enrichment_signTest.Rdata ~/.
  # scp aob2x@rivanna.hpc.virginia.edu:~/drosRtec_enrichment_cline.Rdata ~/.
  load("~/dest_core20.enrichment_signTest.Rdata")
  load("~/drosRtec_enrichment_cline.Rdata")


### basic plot
  en <-
  ggplot(data=o.dest_core20[thr==0.05][mod=="core20_DEST"]) +
  geom_hline(aes(yintercept=1)) +
  geom_linerange(aes(x=chr, ymin=lci, ymax=uci, group=inv, color=inv), position=position_dodge(width = 1)) +
  geom_point(aes(x=chr, y=or, group=inv, color=inv), position=position_dodge(width = 1)) +
  theme_bw() +
  ylab("Odds ratio") +
  ggtitle("Enrichment of SNPs in top 5%\nCville GLM vs. Core20(noVA)")


  concord <-
  ggplot(data=o.dest_core20[thr==0.05][mod=="core20_DEST"]) +
  geom_hline(aes(yintercept=.5)) +
  geom_linerange(aes(x=chr, ymin=st.lci, ymax=st.uci, group=inv, color=inv), position=position_dodge(width = 1)) +
  geom_point(aes(x=chr, y=st.pr, group=inv, color=inv), position=position_dodge(width = 1)) +
  theme_bw() +
  ylab("Concordance Rate") +
  ggtitle("Concordance of allele frequency\nchange of SNPs in top 5%\nCville GLM vs. Core20(noVA)")

  en.cline <-
  ggplot(data=o.cline[thr==0.05][mod=="cline"]) +
  geom_hline(aes(yintercept=1)) +
  geom_linerange(aes(x=chr, ymin=lci, ymax=uci, group=inv, color=inv), position=position_dodge(width = 1)) +
  geom_point(aes(x=chr, y=or, group=inv, color=inv), position=position_dodge(width = 1)) +
  theme_bw() +
  ylab("Odds ratio") +
  ggtitle("Enrichment of SNPs in top 5%\nCville GLM vs. Clinal (Machado)")


  concord.cline <-
  ggplot(data=o.cline[thr==0.05][mod=="cline"]) +
  geom_hline(aes(yintercept=.5)) +
  geom_linerange(aes(x=chr, ymin=st.lci, ymax=st.uci, group=inv, color=inv), position=position_dodge(width = 1)) +
  geom_point(aes(x=chr, y=st.pr, group=inv, color=inv), position=position_dodge(width = 1)) +
  theme_bw() +
  ylab("Concordance Rate") +
  ggtitle("Concordance of allele frequency\nchange of SNPs in top 5%\nCville GLM vs. Clinal (Machado")



layout <- "
AB
CD"


  mega <- en + concord +  en.cline + concord.cline + plot_layout(design=layout, guides = 'collect')

  ggsave(mega, file="Overwintering_18_19/Core20_enrichment/DEST_core20/core20_enrich_concord.pdf", h=8, w=8)
