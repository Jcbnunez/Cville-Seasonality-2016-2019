rm(list = ls())

library(tidyverse)

#### output files:
glm_models_out <-  "/project/berglandlab/alan/dest_glm_final_nested_qb_alternatives_AIC/dest_glm_final_nested_qb_alternatives_AIC"

#system( paste("ls ", glm_models_out, " | head", sep = "") )


load("/project/berglandlab/alan/dest_glm_final_nested_qb_alternatives_AIC/dest_glm_final_nested_qb_alternatives_AIC/job100.Rdata")