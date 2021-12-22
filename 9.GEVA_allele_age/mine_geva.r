#Extract age of alleles


library(tidyverse)
library(data.table)
library(magrittr)

load("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata") 


glm.out %>% 
filter(
	  #rnp.clean < 0.01,
	  mod == "aveTemp+year_factor") %>% 
mutate(SNP_id = paste(chr, pos, sep = "_") ) ->
glm.out


geva_age <- fread("/scratch/yey2sn/Overwintering_ms/9.GEVA/Dmel_GEVA.0.01.TMRCA.txt")

geva_age %<>%
separate(V1, into = c("chr", "array_start", "array_end")) %>%
mutate(SNP_id = paste(chr, V2, sep = "_") )

left_join(glm.out, geva_age) -> glm.out.geva

glm.out.geva %>% 
filter(V8 > 1) %>%
filter(chr == "2L") %>%
ggplot(aes(
y=-log10(rnp.clean),
x=V8
)) +
geom_point() ->
geva_glm

ggsave(geva_glm, file = "geva_glm.png")