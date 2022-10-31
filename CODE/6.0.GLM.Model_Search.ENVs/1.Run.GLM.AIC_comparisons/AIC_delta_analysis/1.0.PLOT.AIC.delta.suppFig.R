### download
#  system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/environmental_ombibus_global/bestAIC.global.Rdata ~/.")

rm(list = ls())

### libraries
library(data.table)
#library(gdata)
library(foreach)
library(doMC)
library(reshape2)
registerDoMC(5)
library(tidyverse)
library(magrittr)

###
load("./nasa_power.weather.mod.Rdata")
names(weather.ave)[1] = "sampleId"
clusts <- fread("./pop.clusters.txt")
full_join(weather.ave, clusts) -> weather.ave.clust

weather.ave.clust$Continental_clusters[grep("VA", weather.ave.clust$sampleId) ] = "5.Cville"

weather.ave.clust %>%
  filter(Continental_clusters %in% c("1.Europe_W", "3.Europe_E", "5.Cville") ) %>%
  filter(!is.na(temp.max)) %>%
  group_by(Continental_clusters) ->
  weather.data.sel

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

### combine with snp identity
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
setkey(snp.dt, id)

#### load the models -->
load("make.Fig2B.dat.Rdata")
#o2.ag %>%
#  filter(!(cluster %in% 
#             c("2.North_America_E", 
#               "2.North_America_I95", 
#               "2.North_America_Mid", 
#               "2.North_America_W")) ) %>%
#  filter(!is.na(sig)) %>%
#  #filter(sig == TRUE) %>%
#  group_by(chr, inv, cluster) %>%
#  arrange(-rr) %>% 
#  mutate(model.rank = 1:n()) %>%
#  mutate(mod.nam = paste(var,mod, cluster, sep = ";")) ->
#  best.nominal.models 
#  
#best.nominal.models %>%
#  filter(model.rank == 1) %>%
#  select(chr, inv, cluster, best.model = mod.nam) ->
#  best.mods
#
#best.nominal.models %>%
#  left_join(best.mods) %>% 
#  filter(mod.nam != best.model) ->
#  best.nominal.models.compared
#
#best.nominal.models.compared %<>%
#  mutate(cor.mods = NA,
#         mod.class = NA,
#         rank.mod.rel.best = NA)
#
#### make weather corrs
#foreach(i = 1:dim(best.nominal.models.compared)[1])%do%{
#  
#  mod1 = best.nominal.models.compared$mod.nam[i]
#  mod2 = best.nominal.models.compared$best.model[i]
#  mod.ran.rel = best.nominal.models.compared$model.rank[i]
#  
#  if(str_split(mod1, ";")[[1]][1] %in% c("pop_year") | 
#     str_split(mod2, ";")[[1]][1] %in% c("pop_year") ){
#    best.nominal.models.compared$cor.mods[i] = NA
#    best.nominal.models.compared$mod.class[i] = "year"
#    best.nominal.models.compared$rank.mod.rel.best[i] = mod.ran.rel
#    
#  } else  if(str_split(mod1, ";")[[1]][1] %in% c("null") | 
#             str_split(mod2, ";")[[1]][1] %in% c("null") ){
#    best.nominal.models.compared$cor.mods[i] = NA
#    best.nominal.models.compared$mod.class[i] = "null"
#    best.nominal.models.compared$rank.mod.rel.best[i] = mod.ran.rel
#  } else{
#    
#    dat1 = weather.data.sel %>%
#      filter(mod == str_split(mod1, ";")[[1]][2]  &
#               Continental_clusters == str_split(mod1, ";")[[1]][3]  ) %>%
#      dplyr::select(str_split(mod1, ";")[[1]][1])
#    names(dat1)[2] = "X"
#    
#    dat2 = weather.data.sel %>%
#      filter(mod == str_split(mod2, ";")[[1]][2]  &
#               Continental_clusters == str_split(mod2, ";")[[1]][3]  ) %>%
#      dplyr::select(str_split(mod2, ";")[[1]][1])
#    names(dat2)[2] = "X"
#    
#    cor.test(x=dat1$X, y=dat2$X)
#    
#    best.nominal.models.compared$cor.mods[i] = cor.test(x=dat1$X, y=dat2$X)$estimate
#    best.nominal.models.compared$mod.class[i] = "env"
#    best.nominal.models.compared$rank.mod.rel.best[i] = mod.ran.rel
#  }
#  
#}
##
##save(best.nominal.models.compared, file = "best.nominal.models.compared.Rdata")
#

best.nominal.models.compared <- get(load("best.nominal.models.compared.Rdata"))

### Prepare model extraction 
root.folder = "/project/berglandlab/alan/environmental_ombibus_global"


#i = 1
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
message(i)
#

               message( paste(i, dim(best.nominal.models.compared)[1], sep = " | " ))      
                #print(best.nominal.models[i,])  
                  
                  file.path_mod1 <- paste(
                    root.folder,
                    best.nominal.models.compared$best.model[i],
                    paste(best.nominal.models.compared$best.model[i], "glmRNP.Rdata", sep = "."),
                    sep = "/"
                  )
            
                  file.path_mod2 <- paste(
                    root.folder,
                    best.nominal.models.compared$mod.nam[i],
                    paste(best.nominal.models.compared$mod.nam[i], "glmRNP.Rdata", sep = "."),
                    sep = "/"
                  )
                
                  #file.path_mod3 <- paste(
                  #  root.folder,
                  #  best.nominal.models$`3`[i],
                  #  paste(best.nominal.models$`3`[i], "glmRNP.Rdata", sep = "."),
                  #  sep = "/"
                  #)
                  
          ### load best mod
          best.mod <- get(load(file.path_mod1))
          best.mod %>%
            dplyr::select(chr,pos,perm, rnp, AIC_best=AIC) ->
            bes.mod
          
          second.mod <- get(load(file.path_mod2))
          second.mod %>%
            dplyr::select(chr,pos,perm, AIC_sec=AIC) ->
            sec.mod
          
          #third.mod <- get(load(file.path_mod3))
          #third.mod %>%
          #  dplyr::select(chr,pos,perm, AIC_thi=AIC) ->
          #  thi.mod
          
          #join
          left_join(bes.mod, sec.mod ) %>%
            #left_join(thi.mod) %>%
            left_join(snp.dt[,c("chr","pos","invName")]) ->
            joint.mods
          
             
          joint.mods %>%
            mutate(aic.1v2 = abs(AIC_best-AIC_sec),
                   #aic.2v3 = abs(AIC_sec-AIC_thi),
                     ) ->
            joint.mods.aic
          
          joint.mods.aic %>%
            filter(!is.na(rnp)) %>%
            mutate(rnp.trsh = case_when(
              rnp >= 0.05 ~ "noRNP",
              rnp < 0.05 & rnp >= 0.01 ~ "rnp5%",
              rnp < 0.01 & rnp >= 0.001 ~ "rnp1%",
              rnp < 0.001  ~ "rnp01%")) %>% 
            mutate(invName.AIC = case_when(
              invName == "none" ~ "Outside Inversion",
              invName != "none" ~ "Inside Inversion",
            )) %>% 
            filter(perm == 0) %>%
            group_by(rnp.trsh, invName.AIC, 
                     chr, perm ) %>% 
            summarise(
                      mean.AICd.1 = mean(aic.1v2),
                      med.AICd.1 = quantile(aic.1v2, 0.5),
                      uci.AICd.1 = quantile(aic.1v2, 0.975),
                      lci.AICd.1 = quantile(aic.1v2, 0.025)
                      #mean.AICd.2 = mean(aic.2v3),
                      ) %>% 
            mutate(
                   cluster = best.nominal.models.compared$cluster[i],
                   chr.focal = best.nominal.models.compared$chr[i],
                   inv.focal = best.nominal.models.compared$inv[i],
                   best.mod = best.nominal.models.compared$best.model[i],
                   other.mod = best.nominal.models.compared$mod.nam[i],
                   cor.mods = best.nominal.models.compared$cor.mods[i],
                   mods.rank.typ = best.nominal.models.compared$mod.class[i],
                   mods.rank.diff = best.nominal.models.compared$rank.mod.rel.best[i],
                   ) ->
            aic.d.ag

         ## aic.d.ag %>%
         ##   filter(invName.AIC == inv.focal) %>%
         ##   filter(chr == chr.focal) ->
         ##   aic.d.ag.flt
          
          ##aic.d.ag.flt %>%
          ##  ggplot(aes(
          ##    x=rnp.trsh,
          ##    y=med.AICd.1,
          ##    ymin=lci.AICd.1,
          ##    ymax=uci.AICd.1,
          ##  )) +
          ##  geom_errorbar() +
          ##  geom_point() +
          ##  facet_wrap(invName.AIC~chr) ->
          ##  aic.d.ag.flt.plot
          ##ggsave(aic.d.ag.flt.plot, file = "aic.d.ag.flt.plot.png")
          
          
file.name = paste("all.aic.test", best.nominal.models.compared$cluster[i], 
                  best.nominal.models.compared$best.model[i],
                  best.nominal.models.compared$mod.nam[i], "Rdata", sep = ".")
          
save(aic.d.ag, file = paste("out.aic.folder", file.name, sep = "/"))
