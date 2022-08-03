### Make Panel B

rm(list = ls())

library(tidyverse)
library(reshape2)
library(magrittr)
library(foreach)
library(viridis)
library(data.table)
library(ggridges)
library(broom)
library(ade4)


load("./FST.geographical.Rdata")

fst.dat.geo.EC$SNP.set %>% unique() -> sets


clusters_to_choose <- c("3.Europe_E", "1.Europe_W")
clusts.mantel.geo = 
    foreach(i=1:length(clusters_to_choose), .combine = "rbind")%do%{
      
      for(k in 1:length(sets) ){
        fst.dat.geo.EC %>%
          filter(SNP.set == sets[k]) %>%
          filter(Continental_clusters == clusters_to_choose[i]) %>%
          dplyr::select(samp1, samp2, FST) -> tmp.raw.top
        
        ##tmp.3 = as.dist(tmp.2)
        tmp.3 = dist(tmp.raw.top$FST)
        #rownames(tmp.3) = tmp.2[1]
        #colnames(tmp.3) = names(tmp.2)[-1]
        assign(paste(sets[k], "matrix", sep = "."), tmp.3)
      }
      
      mantel.rtest(glm.snps.matrix, macthed.controls.noInv.matrix, 
                   nrepet = 999, alt = "two-sided") -> glm.con1
      mantel.rtest(glm.snps.matrix, macthed.controls.Inv.matrix, 
                   nrepet = 999, alt = "two-sided") -> glm.con2
      mantel.rtest(macthed.controls.Inv.matrix, macthed.controls.noInv.matrix, 
                   nrepet = 999, alt = "two-sided") -> conts
      
      data.frame(
        pop = rep(clusters_to_choose[i], 3),
        comp = c("glm-ctrNoInv", "glm-ctrInv", "ctrNoInv-ctrInv"),
        obs = c(glm.con1$obs , glm.con2$obs , conts$obs ),
        p.vals = c(glm.con1$pvalue, glm.con2$pvalue, conts$pvalue)
      )
  }

clusts.mantel.geo %<>% mutate(sig = case_when(p.vals > 0.01 ~ "Different",
                                          p.vals <= 0.01 ~ "Same"))

write.table(clusts.mantel.geo, file = "geo.mantelTest.test.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
