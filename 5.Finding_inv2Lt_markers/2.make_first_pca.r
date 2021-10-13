#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(RVenn)

########
######## Now build PCA
########

load("./imputated_data_DGRP_dat.2L.Rdata")


imputated_data_DGRP[,grep("STD", colnames(imputated_data_DGRP))] %>%
  .[,which(colnames(.) != "STD_line_48" )] -> std_train

imputated_data_DGRP[,grep("INV", colnames(imputated_data_DGRP))]  -> inv_train

#viz pca
cbind(std_train,
      inv_train
      ) -> train_subset

##Save individuals to use

write.table(names(train_subset), 
            file = "./lines_used_for_2lt_demarcation.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)


train_subset %>% 
  t() %>%
  as.data.frame() %>% 
  PCA(scale.unit = FALSE, 
      ncp = 5,
      graph = F
      ) -> 
  PCA_obj_2LT

save(PCA_obj_2LT, file = "PCA_obj_2LT.Rdata")

data.frame(PCA_obj_2LT$ind$coord[,1:2]) %>%
  as.data.frame() %>%
  mutate(ind = row.names(.) ) %>% 
  separate(ind, into = c("invst", "line","lineid") , sep = "_", remove = F) %>%
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.2,
      color = invst
    )
  ) + geom_point(size = 3) ->
  inv_pca_fig

ggsave(inv_pca_fig,
       file = "inv_pca_fig.pdf",
       width = 4,
       height  = 3)

