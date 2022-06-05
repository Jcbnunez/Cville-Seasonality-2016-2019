### libraries
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### load SNP level tests
  load("~/glm.out.lrt.Rdata")


### load weather data
  load(file="~/weather_impute.Rdata")

  wm.impute[,term:=paste(variable, dateDelta, sep="_")]

### samps
  samps <- fread("~/DEST_10Mar2021_POP_metadata.csv")

  samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
  samps[nchar(month)==1, month:=paste("0", month, sep="")]
  samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
  samps[nchar(day)==1, day:=paste("0", day, sep="")]

  setkey(samps, sampleId)
  samps <- samps[J(unique(wm.impute$sampleId))]

### weather averages
  temps.summary <- wm.impute[J(samps[locality=="VA_ch"][year>2012]$sampleId)][variable=="tmax"][,list(mu=mean(value)), list(sampleId)]

################
### SNP data ###
################

  ### open GDS for common SNPs (PoolSNP)
   genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
