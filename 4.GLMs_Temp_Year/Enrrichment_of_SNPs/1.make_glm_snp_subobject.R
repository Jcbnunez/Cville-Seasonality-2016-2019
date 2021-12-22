#Load libraries
library(data.table)
library(tidyverse)
library(reshape2)
library(magrittr)

#Load data
#######################
## load the annotation file
annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
load(annotation_file)

# add link to data --
loaded_data_list = list()

reps = 0:100
for(i in 1:2){
  print(i)
  tmp <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", reps[i], ".Rdata", 
               sep = "")
  print(tmp)
  load(tmp)
  
  left_join(glm.out, 
            annotation_dmel_sp) ->
    glm_data_ALL_annotated
  
  glm_data_ALL_annotated %>%
    separate(Consequence, into = c("Consequence", "etc"),
             sep = "_") -> glm_Cats
  
  loaded_data_list[[i]] = glm_data_ALL_annotated
}

glm_data = do.call(rbind, loaded_data_list)

#glm_data_ALL_annotated$Consequence %>% table

#######################
# Check behaviour

save(glm_data,
     file = "glm_data.annotated.Rdata")  
