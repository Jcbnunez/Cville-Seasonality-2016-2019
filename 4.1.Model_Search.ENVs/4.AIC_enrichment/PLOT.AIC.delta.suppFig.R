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

### combine with snp identity
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")
setkey(snp.dt, id)

#### load the models -->
load("make.Fig2B.dat.Rdata")
o2.ag %>%
  filter(!(cluster %in% 
             c("2.North_America_E", 
               "2.North_America_I95", 
               "2.North_America_Mid", 
               "2.North_America_W")) ) %>%
  filter(!is.na(sig)) %>%
  filter(sig == TRUE) %>%
  group_by(chr, inv, cluster) %>%
  arrange(-rr) %>%
  slice_head(n=3) %>%
  mutate(model.rank = 1:n()) %>%
  as.data.frame() %>%
  mutate(mod.nam = paste(var,mod, cluster, sep = ";")) %>% 
  reshape2::dcast(cluster+chr+inv~model.rank, value.var = "mod.nam")  ->
  best.nominal.models

### Prepare model extraction 
root.folder = "/project/berglandlab/alan/environmental_ombibus_global"

models.ext.ag <- foreach(i= 1:dim(best.nominal.models)[1] , 
                .errorhandling="remove")%dopar%{
                
               message( paste(i, dim(best.nominal.models)[1], sep = " | " ))      
                #print(best.nominal.models[i,])  
                  
                  file.path_mod1 <- paste(
                    root.folder,
                    best.nominal.models$`1`[i],
                    paste(best.nominal.models$`1`[i], "glmRNP.Rdata", sep = "."),
                    sep = "/"
                  )
            
                  file.path_mod2 <- paste(
                    root.folder,
                    best.nominal.models$`2`[i],
                    paste(best.nominal.models$`2`[i], "glmRNP.Rdata", sep = "."),
                    sep = "/"
                  )
                
                  file.path_mod3 <- paste(
                    root.folder,
                    best.nominal.models$`3`[i],
                    paste(best.nominal.models$`3`[i], "glmRNP.Rdata", sep = "."),
                    sep = "/"
                  )
                  
          ### load best mod
          best.mod <- get(load(file.path_mod1))
          best.mod %>%
            dplyr::select(chr,pos,perm, rnp, AIC_best=AIC) ->
            bes.mod
          
          second.mod <- get(load(file.path_mod2))
          second.mod %>%
            dplyr::select(chr,pos,perm, AIC_sec=AIC) ->
            sec.mod
          
          third.mod <- get(load(file.path_mod3))
          third.mod %>%
            dplyr::select(chr,pos,perm, AIC_thi=AIC) ->
            thi.mod
          
          #join
          left_join(bes.mod, sec.mod ) %>%
            left_join(thi.mod) %>%
            left_join(snp.dt[,c("chr","pos","invName")]) ->
            joint.mods
          
             
          joint.mods %>%
            mutate(aic.1v2 = abs(AIC_best-AIC_sec),
                   aic.2v3 = abs(AIC_sec-AIC_thi),
                     ) ->
            joint.mods.aic
          
          joint.mods.aic %>%
            mutate(rnp.trsh = case_when(
              rnp < 1.00 & rnp >= 0.05 ~ "noRNP",
              rnp < 0.05 & rnp >= 0.01 ~ "rnp5%",
              rnp < 0.01 & rnp >= 0.01 ~ "rnp1%",
              rnp < 0.01 & rnp <= 0.001 ~ "rnp01%")) %>% 
            group_by(rnp.trsh, invName=="none", 
                     chr ) %>% 
            summarise(mean.AICd.1 = mean(aic.1v2),
                      mean.AICd.2 = mean(aic.2v3),
                      ) %>% 
            mutate(cluster = best.nominal.models$cluster[i],
                   chr.focal = best.nominal.models$chr[i],
                   inv.focal = best.nominal.models$inv[i]
                   ) ->
            aic.d.ag
          
          aic.d.ag %>%
            filter(chr == chr.focal)
                  
                }

