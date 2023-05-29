### Step 1. Prepare the data

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)
library(doMC)
library(lubridate)

registerDoMC(4)

### User provided parameters
args = commandArgs(trailingOnly=TRUE)

### Import some user given parameters
pop = args[1]
print(paste("now loading", pop))
#pop <- "VA_ch"bn14ntlqbeo#
#
p_tresh = args[2]
#p_tresh = 0.05
print(paste("now processing p <", p_tresh))

#One could, also turn the model into an iterator to compare the T and Y+T models
#for now just Y+T
model = "aveTemp+year_factor"

### add root map
root_path <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out."

### Extratc SNPs
iterations = c(0,1:98,100)

#output lists

for(i in 1:length(iterations)){ ### open i
  
  print(iterations[i])
  
  load(paste(root_path,
             pop,
             "_", 
             iterations[i], 
             ".Rdata", sep = ""))
  
  category = ifelse(iterations[i] == 0, "Obs", "Per")
  print(category)
  
  glm.out %>%
    filter(mod == model) %>%
    filter(chr %in% c("2L")) %>%
    filter(complete.cases(rnp.clean)) %>%
    filter(rnp.clean < p_tresh) %>% 
    dplyr::select(chr, pos) %>%
    filter(
      pos >= 2225744, 
      pos <= 13154180
    ) %>%
    mutate(
      feature = "SNP",
      type = "GLM_outlier_inv") ->
    GLM_outlier
  
inv_markers_id <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)
    
inv_markers_id %>%
  separate(V1, into = c("chr","pos", "feature"), sep = "_") %>%
  mutate(type = "Inv2Lt_marker") ->
  inversion_marker
    
rbind(GLM_outlier,
      inversion_marker) %>%
  .[order(as.numeric(.$pos)),] %>%
  mutate(Variant_id = paste(chr, pos, feature, 
                            sep = "_"),
         Pop = pop,
         category = category,
         ith = iterations[i],
         p_tresh = p_tresh
         ) ->
  inversion_glm_markers


write.table(inversion_glm_markers, 
            file = paste(pop ,p_tresh, paste("ith",iterations[i], sep = ""), "snp_master_list.txt", sep = "_"  )  ,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(inversion_glm_markers$Variant_id, 
            file = paste(pop ,p_tresh, paste("ith",iterations[i], sep = ""), "variant_id.txt", sep = "_"  )  ,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

} ### Close i



