library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
args[2] -> permnum

chropenXTX <- function(chr1){
  pos <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_positiontable.txt"), h=T)
  if (permnum == 0){
  XtXtab <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_locostitch_summary_pi_xtx.out"), h=T)
  }
  else {
    XtXtab <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_perm_", permnum, "_summary_pi_xtx.out"), h=T) 
  }
  dim(XtXtab)[1] -> length
  XtXtab <- data.table("chr" = rep(chr1, length), pos, XtXtab)
  return(XtXtab)
}

chropenBF <- function(chr1){
  pos <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_positiontable.txt"), h=T)
  pos10 <- data.table("pos" = rep(pos$pos, 10), "chr" = rep(pos$rep, 10), "variant.id" = rep(pos$variant.id, 10), "invName" = rep(pos$invName, 10))
  if (permnum == 0){
    BFtab <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_locostitch_summary_betai_reg.out"), h=T)
  } else {
    BFtab <- read.table(paste0("/scratch/nzx3cc/nzx3cc/perms/",chr1,"perms/", chr1, "_perm_", permnum, "_summary_betai_reg.out"), h=T)
  }
  dim(BFtab)[1] -> length
  BFtab <- data.table("chr" = rep(chr1, length), pos10, BFtab)
  return(BFtab)
}
chr <- c("2L", "2R", "3L", "3R")
xtx <- lapply(chr, chropenXTX)
bf <- lapply(chr, chropenBF)

xtx <- rbindlist(lapply(xtx, as.data.frame.list))
bf <- rbindlist(lapply(bf, as.data.frame.list))

setkey(xtx, pos)
setkey(bf, COVARIABLE, pos)

if (permnum == 0){
  write.table(xtx, "wholeLOCOXtX.out", col.names=T, quote=F)
  write.table(bf, "wholeLOCOBF.out", col.names=T, quote=F)
} else {
  write.table(xtx, paste0("wholeLOCO_perm_",permnum, "_XtX.out"), col.names=T, quote=F)
  write.table(bf, paste0("wholeLOCO_perm_",permnum, "_BF.out"), col.names=T, quote=F)
}
