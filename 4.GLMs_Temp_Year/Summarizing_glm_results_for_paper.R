## summarizing SNPs

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

###

library(tidyverse)
library(magrittr)

##
##
## Load data:

load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))

glm.out %>%
  filter(mod == "year_factor") %>%
  filter(p.lrt < 0.01) %>% 
  dim %>% .[1] -> out_alph1

glm.out %>%
  filter(mod == "year_factor") %>% dim %>% .[1] -> out_all

out_alph1/out_all -> perc_real


# Load real data

percent_list= list()
for(i in 1:100){
  
  if(i == 99){break}
  
load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))

glm.out %>%
filter(mod == "year_factor") %>%
filter(p.lrt < 0.01) %>% 
  dim %>% .[1] -> out_alph1

glm.out %>%
filter(mod == "year_factor") %>% dim %>% .[1] -> out_all

out_alph1/out_all -> perc_perm
percent_list[[i]] = perc_perm
} # close i


