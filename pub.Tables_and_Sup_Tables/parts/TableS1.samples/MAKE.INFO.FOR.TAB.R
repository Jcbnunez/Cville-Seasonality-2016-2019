library(vroom)
library(tidyverse)
library(data.table)
library(magrittr)

##
chosen = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/samps.chosen.for.CM.txt")
sras = vroom("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/SRA.alyssa.samps.txt")

sras %>%
  filter(Sample.Name %in%  chosen$Sample.Name ) %>%
  arrange(Sample.Name) ->
  used.samps

write.table(used.samps, file = "/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/used.samps.cm.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


######
######
######

#tab1 <- fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/TableS1.samps.pooled.ind.txt")

samps <- read.csv("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/samps.csv")
samps.used <- fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/samps.used.inpaper.txt")


left_join(samps.used, samps) %>% 
  dplyr::select(sampleId, sampleType) ->
  sample.typeDEST

write.table(sample.typeDEST, 
            file = "/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/sample.typeDEST.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

####
#vcf="CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz"
#module load bcftools
#bcftools query -l $vcf

priscdat = read.csv("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/Priscilla/overwintering_DNA_plates_051920.csv")
used.samps.p = fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/Priscilla/samps.used.txt")

priscdat %<>%
  mutate(sampleId = paste(Plate, Well, paste("CM", Season, sep = ""), sep = "_") )

left_join(used.samps.p, priscdat) %>%
  mutate(Date = as.Date(Date, format = "%m.%d.%Y")) ->
  used.samps.prisd

write.table(used.samps.prisd, 
            file = "/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/pub.Tables_and_Sup_Tables/parts/TableS1.samples/used.samps.prisd.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



