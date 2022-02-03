library(tidyverse)
library(data.table)


#pop = args[1]
#print(paste("now loading", pop))
pop <- "VA_ch"
#p_tresh = args[2]
p_tresh = 0.05
#print(paste("now processing p <", p_tresh))

root_path <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out."

iterations = c(0,1:98,100)
i=1

load(paste(root_path,
           pop,
           "_", 
           iterations[i], 
           ".Rdata", sep = ""))


chrs_selection = list(
                      all=c("2L","2R","3L", "3R"),
                      chr2L="2L",
                      chr2R="2R",
                      chr3L="3L",
                      chr3R="3R")

## Loop over chromosomes
for(k in 1:length(chrs_selection)){
  
  glm.out %>%
    filter(
      mod == "aveTemp+year_factor",
      chr %in% chrs_selection[[k]] )  %>%
    .[complete.cases(.$rnp.clean),] %>%
    dplyr::select(chr, pos) ->
    glm_universe
  
  fwrite(glm_universe, 
         file = paste(names(chrs_selection)[[k]],"noTresh", "glm_universe.txt", sep = "_"), 
         append = FALSE, 
         quote = FALSE, 
         sep = "\t",
         row.names = FALSE,
         col.names = FALSE)

  glm.out %>%
    .[complete.cases(.$rnp.clean),] %>%
    filter(
      rnp.clean<p_tresh,
      mod == "aveTemp+year_factor",
      chr %in% chrs_selection[[k]] )  %>%
    dplyr::select(chr, pos) ->
    glm_outliers
  
  fwrite(glm_outliers, 
         file = paste(names(chrs_selection)[[k]],p_tresh, "glm_outliers.txt", sep = "_"), 
         append = FALSE, 
         quote = FALSE, 
         sep = "\t",
         row.names = FALSE,
         col.names = FALSE)
  
} # close K

#save a text file
#########


