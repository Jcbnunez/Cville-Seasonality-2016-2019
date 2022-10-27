### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load data
  load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/2.InversionEnrichmentAnalysis/inversion_enrichment_output.Rdata")

### summarize
  invEnv.ag <- invEnv[,list(or=(TT[perm==0]/TF[perm==0])/(TT[perm!=0]/TF[perm!=0]), perm=perm[perm!=0]), list(locality, i, factor)]
  invEnv.ag.ag <- invEnv[,list(or=mean(log2(or)), sd=sd(log2(or))), list(locality, i, factor)]



  ggplot() +
  geom_violin(data=invEnv[i==0.01][perm!=0],
          aes(x=factor, y=log2(or))) +
  geom_point(data=invEnv[i==0.01][perm==0],
          aes(x=factor, y=log2(or)), color="red") +
  facet_grid(locality~.)


  ggplot(data=invEnv[i==0.01][factor=="2Lt"][locality=="VA_ch"][perm<=30]) +
  geom_point(aes(x=as.numeric(as.factor(or)), y=log2(or), color=as.factor(perm==0))) +
  geom_linerange(aes(x=as.numeric(as.factor(or)), ymin=log2(or.lci), ymax=log2(or.uci), color=as.factor(perm==0)))



mean(
invEnv[i==0.01][factor=="2Lt"][locality=="VA_ch"][perm==0]$or <
invEnv[i==0.01][factor=="2Lt"][locality=="VA_ch"][perm!=0]$or)



          invEnv[i==0.01][perm==0][locality=="VA_ch"]

### basic plot
  gg.bin <- ggplot(data=invEnv[order(-perm)],
              aes(x=log10(i), y=log2(or), color=as.factor(perm==0), group=perm)) +
  geom_line() +
  facet_grid(factor~locality) +
  theme(axis.text.x=element_text(angle=90)) +
  theme_bw()

  gg.bin





  gg <- ggplot(data=invEnv.ag.ag,
              aes(x=log10(i), y=(or))) +
  geom_ribbon(aes(ymin=or-2*sd, ymax=or+2*sd), fill="black", alpha=.25) +
  geom_line() +
  facet_grid(factor~locality) +
  theme(axis.text.x=element_text(angle=90)) +
  theme_bw()
  ggsave(gg, file="~/inversion_enrichment.png", h=11, w=8.5)



















gg.qb <- ggplot(data=o[order(-perm)],
            aes(x=log10(i), y=log2(or), color=as.factor(perm==0), group=perm)) +
geom_line() +
facet_grid(factor~locality) +
theme(axis.text.x=element_text(angle=90)) +
theme_bw()

mega <- gg.bin + gg.qb

ggsave(mega, file="~/fet_out.png", h=11, w=11)




o.ag2 <- o[,list(or=(TF[perm!=0]/TT[perm!=0])/(TF[perm==0]/TT[perm==0]),
                or2=mean(or[perm==0]),
                n_perm=TF[perm!=0] + TT[perm!=0],
                n_obs=TF[perm==0] + TT[perm==0]),
            list(model, factor, locality)]



  o.ag2.ag <- o.ag2[,list(or=mean((or)),
                  or2=mean((or2))),

              list(model, factor, locality)]



ggplot(data=o.ag2.ag) +
 geom_point(aes(x=locality, y=or)) +
 geom_point(aes(x=locality, y=or2), color="red") +
 facet_grid(model~factor) +
 theme(axis.text.x=element_text(angle=90))
