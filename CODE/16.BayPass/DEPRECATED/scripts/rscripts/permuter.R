
library(data.table)

args <- commandArgs(trailingOnly = T)
outprefix <- as.character(args[2])
seed <- args[3]
set.seed(seed)
file <- read.table("/scratch/nzx3cc/nzx3cc/env_analysis/envdata.txt")
permuter <- function(i){
  return(sample(file[i,]))
}
tab <- lapply(c(1:10), permuter)
tab <- rbindlist(lapply(tab, as.data.frame.list))
write.table(tab, file = paste0(outprefix, "_perm_envtable.txt"), row.names=F, col.names=F, quote=F, eol="\n")
file <- read.table("/scratch/nzx3cc/nzx3cc/env_analysis/contrasttable.txt")
tab2 <- lapply(c(1:3), permuter)
tab2 <- rbindlist(lapply(tab2, as.data.frame.list))
write.table(tab2, file = paste0(outprefix, "_perm_contrasttable.txt"), row.names=F, col.names=F, quote=F, eol="\n")

