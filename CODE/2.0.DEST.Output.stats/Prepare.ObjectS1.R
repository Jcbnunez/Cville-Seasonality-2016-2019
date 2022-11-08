#### Prepare object S1.
#### 
#### Please run this script after complteting all analyses in the paper
#### 
#### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90)) %>%
        mutate(time_window = paste(start, end, sep = "_")) %>%
        dplyr::select(mod_id=mod, time_window)




### Load best model
pheno.data = get(load("./SNP.phenos.Rdata"))
pheno.data$pos = as.integer(pheno.data$pos)

###
glm.out = get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"))

glm.out %>% head
names(glm.out)
### summarize data
glm.out %<>%
  mutate(perm.stat = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm"
                               ))
glm.out %>% 
  filter(perm == 0) ->
  glm.out.0
glm.out %>% 
  filter(perm != 0) ->
  glm.out.perms

glm.out.perms %>%
  group_by(chr, pos) %>%
  summarize(Perm.rnp.0.01.quant = quantile(rnp, 0.01)) ->
  glm.out.perms.summ

#### Merge with Perms/Phenos
glm.out.0 %>% 
  dplyr::select(chr, pos, AIC, variable, mod_id=mod, 
                p_lrt, b_temp, se_temp, cluster,cm_mb, invName, rnp ) %>%
  left_join(glm.out.perms.summ) %>% 
  left_join(sets) %>%
  left_join(pheno.data) ->
  glm.out.0.met.VA

####
glm.out.EUE = get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.ave;9;3.Europe_E/temp.ave;9;3.Europe_E.glmRNP.Rdata"))

glm.out.EUE  %>% head
names(glm.out.EUE )
### summarize data
glm.out.EUE  %<>%
  mutate(perm.stat = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm"
  ))

glm.out.EUE  %>% 
  filter(perm == 0) ->
  glm.out.EUE.0
glm.out.EUE  %>% 
  filter(perm != 0) ->
  glm.out.EUE.perms

glm.out.EUE.perms %>%
  group_by(chr, pos) %>%
  summarize(Perm.rnp.0.01.quant = quantile(rnp, 0.01)) ->
  glm.out.EUE.perms.summ

#### Merge with Perms/Phenos

glm.out.EUE.0 %>% 
  dplyr::select(chr, pos, AIC, variable, mod_id=mod, 
                p_lrt, b_temp, se_temp, cluster,cm_mb, invName, rnp ) %>%
  left_join(glm.out.EUE.perms.summ) %>% 
  left_join(sets) %>%
  left_join(pheno.data) ->
  glm.out.EUE.0.met


##### EUW
glm.out.EUW = get(load("/project/berglandlab/alan/environmental_ombibus_global/humidity.ave;8;1.Europe_W/humidity.ave;8;1.Europe_W.glmRNP.Rdata"))

glm.out.EUW  %>% head
names(glm.out.EUW )
### summarize data
glm.out.EUW  %<>%
  mutate(perm.stat = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm"
  ))

glm.out.EUW  %>% 
  filter(perm == 0) ->
  glm.out.EUW.0
glm.out.EUW  %>% 
  filter(perm != 0) ->
  glm.out.EUW.perms

glm.out.EUW.perms %>%
  group_by(chr, pos) %>%
  summarize(Perm.rnp.0.01.quant = quantile(rnp, 0.01)) ->
  glm.out.EUW.perms.summ

#### Merge with Perms/Phenos

glm.out.EUW.0 %>% 
  dplyr::select(chr, pos, AIC, variable, mod_id=mod, 
                p_lrt, b_temp, se_temp, cluster,cm_mb, invName, rnp ) %>%
  left_join(glm.out.EUW.perms.summ) %>% 
  left_join(sets) %>%
  left_join(pheno.data) ->
  glm.out.EUW.0.met

#### joint Rbind

rbind(glm.out.0.met.VA, glm.out.EUE.0.met, glm.out.EUW.0.met) -> glm.best.models

glm.best.models$variable %>% table
glm.best.models$cluster %>% table

##### 
save(glm.best.models, file = "Best.GLM.models.output.Rdata")

names(glm.best.models)


