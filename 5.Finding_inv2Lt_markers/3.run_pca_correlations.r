#Load libraries
library(tidyverse)
library(FactoMineR)

load("./PCA_obj_2LT.Rdata")

dimdesc(PCA_obj_2LT,
        axes = 1, 
        proba = 1
        ) -> PCA_obj_2LT_dimdesc

save(PCA_obj_2LT_dimdesc,
     file = "PCA_obj_2LT_dimdesc.Rdata")
