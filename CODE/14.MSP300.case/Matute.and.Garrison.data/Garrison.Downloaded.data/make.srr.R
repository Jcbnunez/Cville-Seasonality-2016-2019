### create SRR for mauritania study

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)

sra.dat <- vroom("Mauritania.SraExperimentPackage.txt", na = "")
sra.dat %<>%
  .[!is.na(.$"/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession"),] %>%
  .[,c("/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession", "/EXPERIMENT_PACKAGE/RUN_SET/RUN/@alias" )]

names(sra.dat) = c("SRR", "samp")
sra.dat$samp = gsub("RNAseq ", "RNAseq.", sra.dat$samp)
sra.dat$samp = gsub("large insert library PE data", "LargeInsert", sra.dat$samp)
sra.dat$samp = gsub(" ", "_", sra.dat$samp)

sra.dat %>%
  group_by(samp) %>%
  slice_head() -> sra.dat.slice

write.table(sra.dat.slice[,c("samp", "SRR")],
            file = "mauritania.srr.fetch.txt", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"),
            fileEncoding = "")
