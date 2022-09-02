##### ----> extract the In2Lt markers
##### 
##### 

### load modules
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
library(vroom)


##### Preprare model SNPs
base <- "/project/berglandlab/alan/environmental_ombibus_global"

######
######
### ----> loac Cville object ... this is a constant througt
file.cvile <- paste(base, "temp.max;2;5.Cville", "temp.max;2;5.Cville.glmRNP.Rdata", sep = "/" )

### load the glm object. 
print(file.cvile)
out.glm.cvile <- get(load(file.cvile))


out.glm.cvile %>% 
  filter(perm == 0 ) %>%
  filter(rnp < 0.05) %>%
  filter(chr == "2L") %>% 
  mutate( win.name =case_when(
        pos > 2800000 & pos < 3200000 ~ "win_3.1",
        pos > 4470000 & pos < 4870000 ~ "win_4.7",
        pos > 4920000 & pos < 5320000 ~ "win_5.1",
        pos > 6000000 & pos < 6400000 ~ "win_6.1",
        pos > 6600000 & pos < 7000000 ~ "win_6.8",
        pos > 9400000 & pos < 9800000 ~ "win_9.6",
        TRUE ~ "noWin"
    )) -> out.glm.cvile.annot

out.glm.cvile.annot %>%
  filter(win.name != "noWin") %>% 
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_") ) ->
  out.glm.cvile.annot.win

save(out.glm.cvile.annot.win, file = "Cville.GLM.rnp5.snps.Rdata")




