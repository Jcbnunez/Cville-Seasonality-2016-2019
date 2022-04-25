#### Generate population level labels
#### 
library(tidyverse)    # data manipulation and visualization
library(data.table)
library(magrittr)
library(FactoMineR)

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM.toPredict_INV_STD.new.Apr1.2022.Rdata")


write.table(sample_metadata, 
            file = "sample_metadata_CM_world_SNPs.2l.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")



  
  
  
  