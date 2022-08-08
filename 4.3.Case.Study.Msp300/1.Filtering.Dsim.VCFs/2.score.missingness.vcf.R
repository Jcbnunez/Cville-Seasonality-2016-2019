####
####
#### Filtering individuals 
#### 

library(vcfR)
library(adegenet)
library(tidyverse)
library(vroom)

#####
vcf <- read.vcfR("Dsim_all_2L.maf0.01.recode.vcf.gz", verbose = TRUE)
x <- vcfR2genlight(vcf)
x.gen <- tab( x , NA.method =  "asis")

na.per.ind <- apply(x.gen, 1, function(x) sum(is.na(x)))
tot.snps = dim(x.gen)[2] 

data.frame(miss=na.per.ind, 
           tot=tot.snps) %>%
  mutate(sampleId = rownames(.)) %>%
  mutate(Nper = miss/tot) ->
  dsim.metadata

dsim.metadata %>%
  filter(Nper < 0.1) %>%
  .$sampleId ->
  D.sim.samps.to.keep


x.gen %>% colnames -> snps.all
snps.all[grep("2L_4996892", snps.all)]

x.gen[D.sim.samps.to.keep, "2L_4996892_207704"] %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) ->
  mutation_of_interst.stus
  
write.table(D.sim.samps.to.keep, file = "D.sim.samps.to.keep.txt", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(mutation_of_interst.stus, file = "mutation_of_interst.stus.txt", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


###
