# ijob -A berglandlab -c10 -p standard --mem=80G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  .libPaths(c("/scratch/aob2x/biol4559/packages_temp/", .libPaths())); .libPaths()

  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(ggplot2)

### load in GLM to look at overlap stats
  load("/scratch/aob2x/Gio_new/temp.max;2;5.Cville/temp.max;2;5.Cville.v2.glmRNP.Rdata")
  #load("/Users/alanbergland/bf/temp.max;2;5.Cville/temp.max;2;5.Cville.v2.glmRNP.Rdata")

### load in BF replicates
  bf_raw <- foreach(i=1:5)%dopar%{
    message(i)
    tmp <- fread(paste("/standard/vol186/bergland-lab/Gio/tempmax_wholegenome_", i, "_summary_betai_reg.out", sep=""))
    tmp[,rep:=i]
    return(tmp)
  }

  bf_raw <- rbindlist(bf_raw)
  setnames(bf_raw, "BF(dB)", "bf_db")
  bf_real <- bf_raw[,list(bf_db.mean=mean(bf_db), bf_db.median=median(bf_db), bf_db.var=var(bf_db), bf_test=bf_db[rep==1],
                             eBPis.mean=mean(eBPis), eBPis.median=median(eBPis), eBPis.var=var(eBPis), eBPis_test=eBPis[rep==1]), list(MRK)]


  bf_index <- fread("/standard/vol186/bergland-lab/Gio/wholegenome_positiontable.txt")
  bf_real <- merge(bf_real, bf_index, by.x="MRK", by.y="V1")

### merge BF with GLM data
  setkey(bf_real, chr, pos)
  setkey(glm.out, chr, pos)
  m.bf <- merge(bf_real, glm.out[perm==0])

### save real
  save(m.bf, file="~/bayPass_output/bf_real.Rdata")

### load thresholds
  load("~/bayPass_output/bf_threshold.Rdata")

### generate distribution of odds ratios from GLM permutations
  setkey(bf_real, chr, pos)
  setkey(glm.out, perm)
  bf.sim.thr.ag <- bf.sim.thr[,list(bf_db.median=median(bf_db.median)), list(thr)]

  bf.overlap <- foreach(perm.i=c(0:100))%dopar%{
    # perm.i <- 0
    message(perm.i)
    tmp <- glm.out[J(perm.i)]
    setkey(tmp, chr, pos)
    m <- merge(bf_real, tmp)

    m.ag <- foreach(i=1:dim(bf.sim.thr.ag )[1], .combine="rbind")%do%{
      m[,list(TT=sum(bf_db.mean>=bf.sim.thr.ag[i]$bf_db.median & rnp<=.05),
              TF=sum(bf_db.mean< bf.sim.thr.ag[i]$bf_db.median & rnp<=.05),
              FT=sum(bf_db.mean>=bf.sim.thr.ag[i]$bf_db.median & rnp> .05),
              FF=sum(bf_db.mean< bf.sim.thr.ag[i]$bf_db.median & rnp> .05), bf_thr=bf.sim.thr.ag[i]$thr, perm=perm.i),
        list(chr, inv=invName.x!="none")]
    }

    m.ag[,or:=(TT/TF)/(FT/FF)]
    m.ag[order(or)]
    return(m.ag)
  }
  bf.overlap <- rbindlist(bf.overlap)

  ord <- paste(rep(c("2L", "2R", "3L", "3R"), each=2), rep(c(TRUE, FALSE), 4), sep="_")
  bf.overlap[,chr_class:=factor(paste(chr, inv, sep="_"), levels=ord)]
  bf.overlap[is.na(chr_class), chr_class:="simulation"]

### save overlap object
  save(bf.overlap, file="~/bayPass_output/bf_overlap.Rdata")
