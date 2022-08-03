### FST 
### 
### 
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



###save seas data
load("./FST.seasonal.Rdata")

fst.dat.EC$SNP.set %>% unique() -> sets
fst.dat.EC %<>%
  mutate(delta_y = abs(year1-year2))

#### test with mantel
cities_to_choose <- c("Munich", "Akaa", "Broggingen", "Odesa", "Charlottesville")

cities.mantel = 
  foreach(j=0:2, .combine = "rbind", .errorhandling = "remove")%do%{
  foreach(i=1:length(cities_to_choose), .combine = "rbind")%do%{
  
  if(cities_to_choose[i] == "Akaa" & j == 2) { }
    else{
      
    for(k in 1:length(sets) ){
  fst.dat.EC %>%
  filter(SNP.set == sets[k]) %>%
  filter(delta_y == j) %>%
  filter(pop1 == cities_to_choose[i]) %>%
  dplyr::select(samp1, samp2, FST) -> tmp.vals
      tmp.3 = dist(tmp.vals$FST)
      assign(paste(sets[k], "matrix", sep = "."), tmp.3)
    }

    mantel.rtest(glm.snps.matrix, macthed.controls.noInv.matrix, 
             nrepet = 999, alt = "two-sided") -> glm.con1
    mantel.rtest(glm.snps.matrix, macthed.controls.Inv.matrix, 
             nrepet = 999, alt = "two-sided") -> glm.con2
    mantel.rtest(macthed.controls.Inv.matrix, macthed.controls.noInv.matrix, 
             nrepet = 999, alt = "two-sided") -> conts

data.frame(
           pop = rep(cities_to_choose[i], 3),
           delta_y = j,
           comp = c("glm-ctrNoInv", "glm-ctrInv", "ctrNoInv-ctrInv"),
           obs = c(glm.con1$obs , glm.con2$obs , conts$obs ),
           p.vals = c(glm.con1$pvalue, glm.con2$pvalue, conts$pvalue)
           )
    }#else
}
}

cities.mantel %<>% mutate(sig = case_when(p.vals > 0.01 ~ "Different",
                                         p.vals <= 0.01 ~ "Same"))



clusters_to_choose <- c("3.Europe_E", "1.Europe_W")
clusts.mantel = 
  foreach(j=0:2, .combine = "rbind", .errorhandling = "remove")%do%{
    foreach(i=1:length(clusters_to_choose), .combine = "rbind")%do%{
      
      for(k in 1:length(sets) ){
        fst.dat.EC %>%
          filter(SNP.set == sets[k]) %>%
          filter(delta_y == j) %>%
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
        delta_y = j,
        comp = c("glm-ctrNoInv", "glm-ctrInv", "ctrNoInv-ctrInv"),
        obs = c(glm.con1$obs , glm.con2$obs , conts$obs ),
        p.vals = c(glm.con1$pvalue, glm.con2$pvalue, conts$pvalue)
      )
    }
  }

clusts.mantel %<>% mutate(sig = case_when(p.vals > 0.01 ~ "Different",
                                          p.vals <= 0.01 ~ "Same"))

rbind(cities.mantel, clusts.mantel) -> full.matel.analysis

write.table(full.matel.analysis, file = "seasonal.mantelTest.test.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


full.matel.analysis %>% filter(comp ==  "glm-ctrInv")

####


### Box plots
fst.dat.EC %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  mutate(month_diff = abs(month1-month2) ) %>% 
  ggplot(
    aes(x= as.factor(year_diff),
        y=FST,
        color = SNP.set)
  ) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  geom_boxplot(outlier.size =  0.5, size = 0.6) +
  facet_wrap(~Continental_clusters)->
  fst.inv.plot.boxplot

ggsave(fst.inv.plot.boxplot, file = "fst.inv.plot.boxplot.pdf", h = 2.5, w = 6)

### larger boxplot
fst.dat.EC %>%
  mutate(pop.set = case_when(pop1 %in% c("Akaa","Broggingen", "Odesa", "Charlottesville", "Munich") ~ pop1,
                             TRUE ~ Continental_clusters
  )) %>%
  mutate(pop.set = factor(pop.set, levels = c("Charlottesville", "1.Europe_W", "Broggingen", "Munich", "Akaa", "3.Europe_E", "Odesa" ) ) ) %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  mutate(month_diff = abs(month1-month2) ) %>% 
  ggplot(
    aes(x= as.factor(year_diff),
        y=FST,
        color = SNP.set)
  ) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  geom_boxplot(outlier.size =  0.5, size = 0.6) +
  theme(legend.position = "bottom") +
  facet_grid(~pop.set)->
  fst.inv.plot.boxplot.poplevel

ggsave(fst.inv.plot.boxplot.poplevel, file = "fst.inv.plot.boxplot.poplevel.pdf", h = 3, w = 7)
