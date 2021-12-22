# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])

## jobId=1

### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### load inversion
  inv.dt <- fread("~/InversionsMap_hglft_v6_inv_startStop.txt")
  setnames(inv.dt, "chrom", "chr")
  setnames(inv.dt, "stop", "end")
  setkey(inv.dt, start, end)

### load glm.out
  fl <- list.files("/scratch/aob2x/summarized_dest_glm/glm_out", full.names=T)

  load(fl[jobId])

### iterate through rank thresholds
 fet.invs <- foreach(inv.i=c(inv.dt$invName), .combine="rbind", .errorhandling="remove")%do%{
   # inv.i="2Lt"
   foreach(rnp=expand.grid(sapply(c(1:9), function(x) x*10^c(-4:-1)))[,1], .combine="rbind", .errorhandling="remove")%dopar%{
     print(paste(inv.i, rnp, sep=" / "))
     mod.i="aveTemp"
     tab <- table(sig=glm.out$rnp.clean<=rnp,
                 Inv=grepl(inv.i, glm.out$invName))

     temp.inv <- fisher.test(tab)
     data.table(model="aveTemp", factor=inv.i, i=rnp, locality=glm.out[1]$locality, perm=glm.out[1]$perm,
               or=temp.inv$estimate, or.lci=temp.inv$conf.int[1], or.uci=temp.inv$conf.int[2],
               p=temp.inv$p.value,
               FF=tab[1,1], FT=tab[1,2], TF=tab[2,1], TT=tab[2,2])
   }
 }
 fet.invs

### save
  save(fet.invs, file=paste("/scratch/aob2x/summarized_dest_glm/inversionEnrichment/", glm.out[1]$locality, "_", glm.out[1]$perm, ".Rdata", sep=""))
