###
###
###

library(tidyverse)
library(vroom)
library(magrittr)
library(reshape2)


all.samps <- vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/TableS1.samples/DEST_10Mar2021_POP_metadata.csv")
cov.check <- vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/TableS1.samples/TableS1.samps.forPCA.EffCovFiltered.txt")
fst.comps <- vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/19.FST_inversion_pooled/guide_files/master.comp.list.fst.txt")
weather.dat.pops <- get(load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/4.1.Model_Search.ENVs/4.AIC_enrichment/nasa_power.weather.mod.Rdata"))


### cove check
cov.check %>% .$sampleId %>% unique() -> samps.pass.cov
cov.check %>% .$locality %>% unique() -> localities.used

### FST set
fst.comps %>% .$samp1 %>% unique() -> fst.set.1
fst.comps %>% .$samp2 %>% unique() -> fst.set.2
unique(fst.set.1, fst.set.2) -> fst.unique.set

### weather dat set.
weather.dat.pops$sampleId %>% unique() -> glm.pops

### create final set
all.samps %>%
  filter(sampleId %in% unique(c(samps.pass.cov, fst.unique.set, glm.pops)) | locality %in% localities.used) %>% 
  dplyr::select(sampleId, country, city, locality, collectionDate, nFlies, SRA_accession, SRA_experiment, type ) ->
  paper.samps

paper.samps %>%
  mutate(PCA.set = sampleId %in% samps.pass.cov,
         FST.set = sampleId %in% fst.unique.set,
         GLM.set = sampleId %in% glm.pops
         ) ->
  paper.samps.annot


write.table(paper.samps.annot, 
            file = "/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/TableS1.samples/samps.pooled.txt"
            , append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

