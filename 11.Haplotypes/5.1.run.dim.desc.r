library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)

load("./pca_fig_std.DGRP.CM.Rdata")

dimdesc_object = dimdesc(pca_fig_std, proba = 1.0)

save(dimdesc_object, file = "dgrp_cm_dimdesc_object.Rdata")

