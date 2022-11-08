#### Prepare object S1.
#### 
#### Please run this script after complteting all analyses in the paper
#### 
#### 

library(tidyverse)
library(magrittr)
library(reshape2)

### Load best model
glm.out = get(load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"))
pheno.data = get(load("./SNP.phenos.Rdata"))

###
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

#### Merge with Perms
glm.out.0 %>% 
  left_join(glm.out.perms.summ) ->
  glm.out.0.met

##### Merge with Phenos
pheno.data$pos = as.integer(pheno.data$pos)

glm.out.0.met %>%
  left_join(pheno.data) ->
  glm.out.0.met.pheno

save(glm.out.0.met.pheno, file = "VA.GLM.best.model.met.pheno.Rdata")

names(glm.out.0.met.pheno)


