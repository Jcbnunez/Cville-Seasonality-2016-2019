library(tidyverse)
library(data.table)

load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

glm.out %>%
  filter(
    mod == "aveTemp",
    chr == "2L")  %>%
    .[complete.cases(.$rnp.clean),] %>%
  select(chr, pos) ->
  glm_universe
 
#save a text file
fwrite(glm_universe, 
       file = "glm_universe.2L.all.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)
 
 
glm.out %>%
  filter(
    rnp.clean<0.01,
    mod == "aveTemp",
    chr == "2L")  %>%
  select(chr, pos) ->
  glm_p1

fwrite(glm_p1, 
       file = "glm_targets.2L.all.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)

#########
#########
#########

glm.out %>%
  filter(
    mod == "aveTemp",
    chr == "2L",
    invName == "2Lt")  %>%
  .[complete.cases(.$rnp.clean),] %>%
  select(chr, pos) ->
  glm_universe_2lt

#save a text file
fwrite(glm_universe_2lt, 
       file = "glm_universe.2L.inversion.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)


glm.out %>%
  filter(
    rnp.clean<0.001,
    mod == "aveTemp",
    chr == "2L",
    invName == "2Lt")  %>%
  select(chr, pos) ->
  glm_p01_2lt

fwrite(glm_p01_2lt, 
       file = "glm_p01_targets.p001.2L.inversion.txt", 
       append = FALSE, 
       quote = FALSE, 
       sep = "\t",
       row.names = FALSE,
       col.names = FALSE)




