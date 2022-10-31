# ijob -A berglandlab_standard -c5 -p standard --mem=8G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)



fl <- list.files("/project/berglandlab/alan/core20_enrichment/", full.names=T)

o.dest_core20 <- foreach(fl.i=fl, .combine="rbind")%dopar%{
  load(fl.i)
  message(fl.i)
  return(o.dest_core20)
}

save(o.dest_core20, file="~/perm.dest_core20.Rdata")









#
#  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
#  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]
#
#
#  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
#  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]
