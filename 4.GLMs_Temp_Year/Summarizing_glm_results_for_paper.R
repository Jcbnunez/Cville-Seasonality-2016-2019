## summarizing SNPs

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

###

library(tidyverse)
library(magrittr)

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


############ ########## ############
############ 
############ 






