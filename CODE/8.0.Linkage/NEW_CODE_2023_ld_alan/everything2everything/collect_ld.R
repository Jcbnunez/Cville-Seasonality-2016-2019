# ijob -A berglandlab_standard -c1 -p dev --mem=8G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  #library(doMC)
  #registerDoMC(4)

### get files
  fl <- list.files("/scratch/aob2x/everythingLD/bits/", full.names=T)
  fl <- fl[grepl("_", fl)]

  obs <- as.numeric(gsub("job", "", tstrsplit(last(tstrsplit(fl, "/")), "_")[[1]]))
  paste(c(1:999)[!c(1:999)%in%obs], collapse=",")

### load
  o <- foreach(i=fl)%dopar%{
    # i <- fl[1]
    message(i)
    load(i)
    return(o)
  }
  o <- rbindlist(o, fill=T)

### save
  save(o, file="~/pairwise_ld_window_50000_10000.Rdata")
