# curl https://cdn.elifesciences.org/articles/67577/elife-67577-supp1-v2.xlsx \
# --output /Users/alanbergland/Documents/GitHub/Overwintering_18_19/Core20_enrichment/DEST_core20/elife-67577-supp1-v2.xlsx
#

### libraries
  library(data.table)
  library(gdata)
  library(lubridate)

### drosRTEC sample data
  drosRTEC <- read.xls("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/Core20_enrichment/DEST_core20/elife-67577-supp1-v2.xlsx")
  drosRTEC <- as.data.table(drosRTEC)

  drosRTEC[City=="State College"]
  samps[city=="State College"]

  drosRTEC[,Sample:=gsub("PA_sc", "PA_st", Sample)]

### DEST sample data
  samps <- fread("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/DEST_10Mar2021_POP_metadata.csv")
  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]
  samps <- samps[set!="dgn"]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

  samps[locality=="UA_od", locality:="UA_Ode"]

### merge
  samps.drosRTEC <- merge(samps, drosRTEC[Core20=="yes"], by.x="sampleId", by.y="Sample", all.y=T)
  samps.use <- samps.drosRTEC
  samps.use[Population%in%c("KA_to_14", "MA_la_12", "MI_bh_14", "CA_es_12"),concord:=F]
  samps.use[is.na(concord), concord:=T]

  samps.use[concord==T & !grepl("VA", sampleId)]

### save
  save(samps.use, file="Overwintering_18_19/Core20_enrichment/DEST_core20/core20_samps.Rdata")
