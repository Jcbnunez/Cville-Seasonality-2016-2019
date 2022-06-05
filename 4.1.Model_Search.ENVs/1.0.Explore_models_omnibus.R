### Collect and plot 
### 

rm(list = ls())

library(tidyverse)
library(magrittr)
library(data.table)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(corrplot)
#library(gdata)
library(foreach)
library(doMC)
registerDoMC(5)

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


####
wd_core = "/project/berglandlab/alan/environmental_ombibus/"

### Part 1 -- build the super AIC plot
load( paste(wd_core, "bestAIC.Rdata", sep = "/")) 
### p_lrt ==> lrt between model vs year
### min vs year aic.

o.ag %>%
  left_join(sets) %>% 
  mutate(mod_details = paste(var,start,end, sep = "_")) ->
  o.ag.meta

load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")

setnames(snp.dt, "id", "variant.id")
setkey(snp.dt, variant.id)
setkey(o.ag.meta, variant.id)

o.ag.meta <- merge(o.ag.meta, snp.dt[,c("chr", "pos", "variant.id", "invName")])

### summarize
o.ag.ag <- o.ag.meta[,
                list(N=length(variant.id)),
                list(perm, mod, var 
                     ,chr, inv=invName!="none"
                     )]


o2 <- o.ag.ag[,list(N, prop=N/sum(N), totalN=sum(N), mod, var),
              list(perm 
                   ,chr, inv
                   )]

o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                  prop.perm.mu=mean(prop[perm!=0]),
                  prop.perm.lci=quantile(prop[perm!=0], .001),
                  prop.perm.uci=quantile(prop[perm!=0], .999),
                  prop.perm.med=median(prop[perm!=0]),
                  prop.rr=mean(log2(prop[perm==0]/prop[perm!=0])),
                  prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
            list(
                 chr, inv, 
                 mod, var)]
o2.ag[,rr:=prop.real/prop.perm.mu]

o2.ag[order(-prop.real)][]
o2.ag[prop.real>(prop.perm.uci)][order(rr)]
o2.ag[prop.rr-2*prop.sd>0]
o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]


### a bit more manipulation
o2.ag[,sig.aic:=log2(prop.real/prop.perm.uci)>0 | log2(prop.real/prop.perm.lci)<0]

o2.rank <- o2.ag[chr=="2L"][inv==T]
o2.rank[,modRank:=rank(rr)]
setkey(o2.rank, mod, var)
setkey(o2.ag, mod, var)
o2.ag <- merge(o2.ag, o2.rank[,c("mod", "var", "modRank")])
o2.ag[,inv:=ifelse(inv, "Inside Inversion", "Outside Inversion")]
o2.ag[,inv:=factor(inv, levels=c("Inside Inversion", "Outside Inversion"))]

o2.ag %<>%
  left_join(sets) %>% 
  mutate(mod_details = paste(var,start,end, sep = "_")) 


o2.ag %>%
  filter(sig.aic == TRUE, rr > 1 & var != "year" ) %>%
  arrange(chr, inv) %>%
  dplyr::select(chr, inv,  mod_details)

o2.ag %>%
  filter(sig.aic == TRUE, rr > 1)  %>% .$mod_details %>% unique() -> best_models

### plot
#o2.ag %>%
#  left_join(o.rnp.ag.meta) %>%
#  as.data.table() ->
#  o2.ag.rnp.meta

o2.ag %>%
  group_by(paste(chr, inv, sep = "_")) %>%
  arrange(rr) %>%
  mutate(rank_chr = 1:n() ) %>% data.table() -> data.for.plot

data.for.plot %>%
  ggplot() + 
  geom_errorbar(data = filter(data.for.plot, sig.aic == FALSE ),
    aes(
        x=rank_chr,
        y=rr,
        ymin=(prop.real/prop.perm.uci), 
        ymax=(prop.real/prop.perm.lci),
        color = sig.aic), 
    alpha = 0.5
  ) +
  geom_point(data = filter(data.for.plot, sig.aic == FALSE ),
    aes( 
        x=rank_chr,
        y=rr,
        color = sig.aic,
        #shape = compound_sig
    ), alpha = 0.5
  ) +
  geom_errorbar(data = filter(data.for.plot, sig.aic == TRUE ),
                aes(
                  x=rank_chr,
                  y=rr,
                  ymin=(prop.real/prop.perm.uci), 
                  ymax=(prop.real/prop.perm.lci),
                  color = sig.aic), 
                alpha = 1
  ) +
  geom_point(data = filter(data.for.plot, sig.aic == TRUE ),
             aes( 
               x=rank_chr,
               y=rr,
               fill = sig.aic,
               #shape = compound_sig
             ), alpha = 1, color = "black", shape = 21
  ) +
  geom_errorbar(    data = filter(data.for.plot, mod_details == "year_NA_NA"), 
                aes(
                  x=rank_chr,
                  y=rr,
                  ymin=(prop.real/prop.perm.uci), 
                  ymax=(prop.real/prop.perm.lci),
                  color = sig.aic), 
                alpha = 1, color = "purple"
  ) +
  geom_point(
    data = filter(data.for.plot, mod_details == "year_NA_NA"), 
    aes(x=rank_chr,
        y=rr,
        #shape = compound_sig
    ), alpha = 1, shape = 23, color = "black", fill = "purple", size = 2.3
  ) +
  geom_errorbar(    data = filter(data.for.plot, mod_details == "null_NA_NA"), 
                    aes(
                      x=rank_chr,
                      y=rr,
                      ymin=(prop.real/prop.perm.uci), 
                      ymax=(prop.real/prop.perm.lci),
                      color = sig.aic), 
                    alpha = 1, color = "red"
  ) +
  geom_point(
    data = filter(data.for.plot, mod_details == "null_NA_NA"), 
    aes(x=rank_chr,
        y=rr,
        #shape = compound_sig
    ), alpha = 1, shape = 22, color = "black", fill = "red", size = 2.3
  ) +
  scale_color_manual(values = c("grey", "steelblue")) +
  scale_fill_manual(values = c("steelblue")) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  geom_hline(aes(yintercept=1)) + xlab("Model") + 
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(trans = "log2") +
  facet_grid(inv~chr)  -> fig.aic

ggsave(fig.aic, file = "fig.aic.pdf", h = 3, w = 6)


##aic.en.plot <-
##  ggplot(data=plot.in[!is.na(sig)][order(!sig)], aes(color=sig)) +
##  geom_hline(aes(yintercept=1)) + xlab("Model") + theme(axis.text.x=element_blank()) +
##  geom_linerange(aes(x=rank_chr, ymin=(prop.real/prop.perm.uci), ymax=(prop.real/prop.perm.lci)),
##                 position=position_dodge2(width=.5), alpha = 0.5) +
##  geom_point(aes(x=rank_chr, y=(rr), shape=inv), position=position_dodge2(width=.5), size=2) +
##  geom_point(data=plot.in[sig == TRUE], aes(x=rank_chr, y=(rr)), position=position_dodge2(width=.5), shape = 21,  size=2, fill = "green") +
##  facet_grid(inv~chr) +
##  theme_bw() +
##  scale_y_continuous(trans = "log2") +
##  geom_text_repel(data=plot.in[!is.na(sig)][sig==T],
##                  aes(x=rank_chr, y=(rr), label=mod_details),
##                  xlim=c(1,50),
##                  size = 3,
##                  point.padding = 1,
##                  min.segment.length = 0,
##                  max.time = 1, max.iter = 1e5,
##                  box.padding = 1)
##
##ggsave(aic.en.plot, file = "aic.en.plot.pdf", h = 4, w = 7)
##
#### add the RNVP analysis
## o.rnp.ag is the output of rnp analys. loop through files and pulls out results of one model.
## the it rank norms p across genome and applies thresholds and ask. whether there are more snps above threshold
## than expected by perm.
### guide file
load("/project/berglandlab/alan/environmental_ombibus/mod_var.Rdata")

####
o.rnp <- foreach(mvi=mod_var, .combine="rbind")%dopar%{
  message(mvi)
  #mvi <- mod_var[67]
  outDir <- paste("/project/berglandlab/alan/environmental_ombibus/", mvi, sep="")
  
  load(file=paste(outDir, "/", mvi, ".rnp_summary.Rdata", sep=""))
  
  o.rnp.ag[,mod_var:=mvi]
  o.rnp.ag
}


o.rnp.ag <- o.rnp[,list(
                        #rr.ave=mean(log2((pr[perm!=0])/(thr))),
                        #rr.sd=sd(log2((pr[perm!=0])/(thr))),
                        #rr=log2(pr[perm==0]/thr)),
                        #totalN=totalN[perm==0],
                        pr.real=pr[perm==0], 
                        pr.perm.mu=mean(pr[perm!=0]),
                        pr.perm.lci=quantile(pr[perm!=0], .02),
                        pr.perm.uci=quantile(pr[perm!=0], .98),
                        pr.real.mu=mean(pr[perm==0]),
                        pr.perm.med=median(pr[perm!=0]),
                        pr.rr=mean(log2(pr[perm==0]/pr[perm!=0])),
                        pr.sd=sd(log2(pr[perm==0]/pr[perm!=0]))),    
                  list(chr, inv, thr, mod_var)]


o.rnp.ag[,sig.pr:=log2(pr.real/pr.perm.uci)>0 | log2(pr.real/pr.perm.lci)<0]

#o.rnp.ag %>%
#  filter(chr == "2L" & inv == TRUE & mod_var == "temp.ave_4") -> temp4
#sanity check
#temp4[,sig.pr:=log2(pr.real/pr.perm.uci)>0 | log2(pr.real/pr.perm.lci)<0]
#temp4[thr == 0.02]
#temp4[thr == 0.01]

o.rnp.ag %<>%
  mutate(inv = case_when(inv == FALSE ~ "Outside Inversion",
                         inv == TRUE ~ "Inside Inversion")) %>%
  separate(mod_var, remove = F, into = c("var", "mod"), sep = "_") %>%
  mutate(mod = as.numeric(mod)) %>%
  left_join(sets) #-> o.rnp.ag.meta

#o.rnp.ag.meta %>%
#  filter(sig.pr == TRUE & chr == '2L')


######
##load( paste(wd_core, "mov_var.rnp.ag.Rdata", sep = "/")) 
##
##o.rnp.ag %>%
##  separate(mod_var,
##           into = c("var", "mod"),sep = "_", remove = F) %>%
##  mutate(mod= as.numeric(mod)) %>%
##  mutate(test_sig = case_when(rr.ave + 2*rr.sd > 0 & rr.ave - 2*rr.sd > 0 ~ "sig",
##                              rr.ave + 2*rr.sd < 0 & rr.ave - 2*rr.sd < 0 ~ "sig",
##                              TRUE ~ "not.sig"
##                              )) %>%
##  mutate(test_dir = case_when(rr.ave + 2*rr.sd > 0 & rr.ave - 2*rr.sd > 0 ~ "over",
##                              rr.ave + 2*rr.sd < 0 & rr.ave - 2*rr.sd < 0 ~ "under",
##                              TRUE ~ "not.enrriched"
##  )) %>%
##  mutate(inv = case_when(inv == FALSE ~ "Outside Inversion",
##                         inv == TRUE ~ "Inside Inversion")) %>%
##  left_join(sets) -> o.rnp.ag.meta
##
##
##

plot.dat = filter(o.rnp.ag, thr == 0.05, ) %>%
            mutate(model.var = paste(var, start , end, sep = "_")) %>%
            filter(model.var %in% best_models) %>%
            group_by(chr) %>%
            arrange(pr.real) %>%
            mutate(rank.mod = 1:n() )

plot.dat %>% 
  filter(sig.pr == TRUE & pr.real/pr.perm.mu > 1 & inv == "Inside Inversion") %>%
  arrange(pr.real/pr.perm.mu) %>%
  mutate(val = pr.real/pr.perm.mu) %>%
  as.data.frame() -> final_models

write.table(final_models, 
            file = "final_models.txt", append = FALSE, quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



plot.dat %>% filter(model.var == "year_NA_NA")

plot.dat %>%
  ggplot() +
  geom_errorbar(
    aes(x = rank.mod,
        ymin = pr.real/pr.perm.lci,
        ymax = pr.real/pr.perm.uci,
        color = inv),
    position=position_dodge(width=0.6), width = 0.1, size = 0.5
  ) +   
  geom_point(
    aes(x = rank.mod,
        y= pr.real/pr.perm.mu, 
        fill = inv), color = "black",
    position=position_dodge(width=0.6),
    shape = 23, size = 1.9
  ) +   
  geom_point( data = filter(plot.dat, model.var == "year_NA_NA"),
    aes(x = rank.mod,
        y= pr.real/pr.perm.mu, 
        fill = inv), color = "black",
    position=position_dodge(width=0.6),
    shape = 24, size = 1.9
  ) +   
  scale_y_continuous(trans = "log2") +
  #coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(.~chr, ncol =2, scales = "free_x")  -> fig.rr.rnpv

ggsave(fig.rr.rnpv, file = "fig.rr.rnpv.pdf", h = 4, w = 6)

