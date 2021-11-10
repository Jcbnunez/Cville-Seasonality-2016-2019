#Load libraries
library(data.table)
library(tidyverse)
library(reshape2)
library(magrittr)

#Load data

# add link to data --
loaded_data_list = list()

reps = 0:10
for(i in 1:11){
  tmp <- paste("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_", reps[i], ".Rdata", 
               sep = "")
  print(tmp)
  load(tmp)
  
  loaded_data_list[[i]] = glm.out
  
}

glm_data = do.call(rbind, loaded_data_list)

#######################
## load the annotation file

annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
load(annotation_file)

#############
##############
# Add annotation to the whole set
left_join(glm_data, 
          annotation_dmel_sp) ->
  glm_data_ALL_annotated

save(glm_data_ALL_annotated,
     file = "glm_data_ALL_annotated.Rdata")

#######################
load("./glm_data_ALL_annotated.Rdata")

glm_data_ALL_annotated %>% head

glm_data_ALL_annotated$Consequence %>% table

#######################
# Check behaviour

glm_data_ALL_annotated %>%
  filter(perm == 0) %>%
  group_by(Consequence) %>%
  summarize(Min_p = min(rnp.clean, na.rm = T)) ->
  min_rnps

#min_rnps %>% 
#  filter(Min_p < 0.01) %>%
#  .$Consequence %>%
#  .[complete.cases(.)]->
#  selected_cats

glm_data_ALL_annotated %>%
  separate(Consequence, into = c("Consequence", "etc"),
           sep = "_") -> glm_Cats

save(glm_Cats,
     file = "glm_Cats.Rdata")  
