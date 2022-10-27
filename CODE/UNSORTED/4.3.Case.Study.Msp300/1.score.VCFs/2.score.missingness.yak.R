
library(vcfR)
library(adegenet)
library(tidyverse)
library(vroom)

#####
vcf <- read.vcfR("/scratch/yey2sn/Overwintering_ms/msp300.case/out.vcfs/NC_052527.2.genotyped.raw.vcf.gz", verbose = TRUE)
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
snps.all[grep("NC_052527.2_9456177", snps.all)]

#x.gen[samps.to.keep, "2L_4996892_207704"] %>%
#  data.frame(loci=.) %>%
#  mutate(sampleId = rownames(.)) ->
#  mutation_of_interst.stus
#
#mutation_of_interst.stus %>%
#  group_by(loci) %>%
#  summarize(N = n())
#
#write.table(mutation_of_interst.stus, file = "yak.mutation_of_interst.stus", 
#            append = FALSE, quote = F, sep = " ",
#            eol = "\n", na = "NA", dec = ".", row.names = F,
#            col.names = TRUE, qmethod = c("escape", "double"),
#            fileEncoding = "")
#
metadata %>%
  filter(sampleId %in% samps.to.keep) %>%
  slice_head(n=9) %>%
  .$sampleId -> samps.for.tree

write.table(samps.for.tree, file = "samps.for.tree.yak.txt", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

###

