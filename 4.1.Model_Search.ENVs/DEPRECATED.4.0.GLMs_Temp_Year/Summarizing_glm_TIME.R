
###

library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions

## set some generalities
p_tresh=0.05

##
##

######## ---> Year MODEL
######## 
### This line finsd the number of LRT p values of the model for real data
## Load data:
i=0
model="year_factor"

load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))

glm.out %>%
  filter(mod == model) %>%
  filter(p.lrt < p_tresh) %>% 
  dim %>% .[1] -> out_alph1

glm.out %>%
  filter(mod == model) %>% dim %>% .[1] -> out_all

out_alph1/out_all -> perc_real


### This line finsd the number of LRT p values of the model for the permutated data

percent_list= list()
for(i in 1:100){
  
  if(i == 99){break}
  
  load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))
  
  glm.out %>%
    filter(mod == model) %>%
    filter(p.lrt < p_tresh) %>% 
    dim %>% .[1] -> out_alph1
  
  glm.out %>%
    filter(mod == model) %>% dim %>% .[1] -> out_all
  
  out_alph1/out_all -> perc_perm
  percent_list[[i]] = perc_perm
} # close i

print(paste("the mean % of time SNPS in perms is",
            mean(unlist(percent_list)*100)))
print(paste("the SD % of time SNPS in perms is",
            sd(unlist(percent_list)*100)))

