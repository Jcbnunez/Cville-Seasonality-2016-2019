### Prepare Fst sites for analysis
### 
### 

rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(vroom)

### Load master SNP file
load("/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata")

head(snp.dt)

#### load model output
base <- "/project/berglandlab/alan/environmental_ombibus_global"
models = c("temp.max;2;5.Cville")
k=1
file <- paste(base, models[k], paste(models[k],"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)

message(models[k])
out.glm <- get(load(file))

inv.markers = vroom("in2lt_ld_47snps_informative_markers.txt", delim = "\t", col_names = F)
inv.markers %<>% 
  mutate(type = "glm.mrks") %>%
  separate(X1, remove = F, into = c("chr", "pos", "snp"), sep = "_") %>%
  mutate(pos = as.numeric(pos))
inv.markers %>%
  filter(pos > 5000000) %>%
  .$pos %>% mean -> mean.pt.right

inv.markers %>%
  filter(pos < 5000000) %>%
  .$pos %>% mean -> mean.pt.left

inv.markers %<>%
mutate(d.from.left = pos - mean.pt.left,
       d.from.right = pos - mean.pt.right)
       
  
####
####

out.glm %>% 
  filter(perm == 0) %>%
  filter(chr == "2L") %>%
  filter(rnp < 0.05) %>%
  mutate(X1 = paste(chr, pos, "SNP", sep = "_")) %>%
  mutate(type = "outliers") %>%
  mutate(d.from.left = pos - mean.pt.left,
         d.from.right = pos - mean.pt.right) -> 
  outliers

out.glm %>%
  filter(perm == 0) %>%
  filter(chr == "2L") %>%
  filter(rnp > 0.8) %>%
  mutate(X1 = paste(chr, pos, "SNP", sep = "_"))%>%
  mutate(type = "not.outliers") %>%
  mutate(d.from.left = pos - mean.pt.left,
         d.from.right = pos - mean.pt.right) ->  
  not.outliers

### make final table
rbind(
inv.markers[,c("X1", "chr", "pos",  "type","d.from.left", "d.from.right")],
outliers[,c("X1", "chr", "pos",  "type","d.from.left", "d.from.right")],
not.outliers[,c("X1", "chr", "pos",  "type","d.from.left", "d.from.right")]
) -> SNPs.for.fst.analysis

write.table(SNPs.for.fst.analysis[,c("chr","pos")], 
            file = "SNPs.for.fst.analysis.ids.only.txt",
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(SNPs.for.fst.analysis, 
            file = "SNPs.for.fst.analysis.Metadat.txt",
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")
