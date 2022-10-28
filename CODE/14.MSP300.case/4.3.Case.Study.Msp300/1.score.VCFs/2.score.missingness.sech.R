####
####
#### Filtering individuals 
#### 

library(vcfR)
library(adegenet)
library(tidyverse)
library(vroom)

#####
vcf <- read.vcfR("/scratch/yey2sn/Overwintering_ms/msp300.case/out.vcfs/NC_045949.1.genotyped.raw.vcf.gz", verbose = TRUE)
x <- vcfR2genlight(vcf)
x.gen <- tab( x , NA.method =  "asis")

na.per.ind <- apply(x.gen, 1, function(x) sum(is.na(x)))
tot.snps = dim(x.gen)[2] 

data.frame(miss=na.per.ind, 
           tot=tot.snps) %>%
  mutate(sampleId = rownames(.)) %>%
  mutate(Nper = miss/tot) ->
  metadata

metadata %>%
  filter(Nper < 0.1) %>%
  .$sampleId ->
  samps.to.keep


x.gen %>% colnames -> snps.all
snps.all[grep("NC_045949.1_5099990", snps.all)]

x.gen[samps.to.keep, "NC_045949.1_5099990"] %>%
  data.frame(loci=.) %>%
  mutate(sampleId = rownames(.)) ->
  mutation_of_interst.stus
  
mutation_of_interst.stus %>%
  group_by(loci) %>%
  summarize(N = n())

write.table(mutation_of_interst.stus, file = "sech.mutation_of_interst.stus", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

mutation_of_interst.stus %>%
  filter(loci %in% c(0,2)) %>%
  group_by(loci) %>%
  slice_head(n=9) %>%
  .$sampleId -> samps.for.tree
  
write.table(samps.for.tree, file = "samps.for.tree.sech.txt", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

###
