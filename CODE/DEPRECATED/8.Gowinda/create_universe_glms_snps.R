library(tidyverse)
library(data.table)


#pop = args[1]
#print(paste("now loading", pop))
pop <- "VA_ch"
#p_tresh = args[2]
p_tresh = 0.05
#print(paste("now processing p <", p_tresh))
annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
load(annotation_file)

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


###Special Cases
load("./hits_perm_binomial.Rdata")

p_tresh = 0.05

glm.out %>%
  .[complete.cases(.$rnp.clean),] %>%
  left_join(annotation_dmel_sp) %>% 
  filter(
    Symbol %in% hit_at_10,
    rnp.clean<p_tresh,
    mod == "aveTemp+year_factor",
    chr == "2L" )  %>%
  dplyr::select(chr, pos) ->
  glm_2L_outliers_at_10

glm.out %>%
  .[complete.cases(.$rnp.clean),] %>%
  left_join(annotation_dmel_sp) %>% 
  filter(
    Symbol %in% hit_at_5,
    rnp.clean<p_tresh,
    mod == "aveTemp+year_factor",
    chr == "2L" )  %>%
  dplyr::select(chr, pos) ->
  glm_2L_outliers_at_5

fwrite(glm_2L_outliers_at_10, 
       file = "glm_2L_outliers_at_10.glm_outliers.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)

fwrite(glm_2L_outliers_at_5, 
       file = "glm_2L_outliers_at_5.glm_outliers.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)

#save a text file
#########


