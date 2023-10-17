# ijob -A berglandlab_standard -c5 -p standard --mem=5G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(5)

### get files
  mainDir <- "/scratch/aob2x/environmental_ombibus_global"
  outDir <- paste(mainDir, "crossClusterEnrichment", sep="/")

  fl <-  list.files(outDir, full.names=T, recursive=T)

  m.ag <- foreach(fl.i=fl)%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    return(m.ag)
  }
  m.ag <- rbindlist(m.ag)

### summarize
  m.ag.ag <- m.ag[,list(crossCluster.or.pr=mean(or[perm==0]>or[perm!=0], na.rm=T),
                        crossCluster.or=or[perm==0],
                        crossCluster.st.pr=mean(st.pr[perm==0]>st.pr[perm!=0], na.rm=T),
                        crossCluster.st=st.pr[perm==0],
                        focalCluster.or.pr=mean(focalT[perm==0]/(focalT[perm==0]+focalF[perm==0]) > focalT[perm!=0]/(focalT[perm!=0]+focalF[perm!=0])),
                        testCluster.or.pr=mean(testerT[perm==0]/(testerT[perm==0]+testerF[perm==0]) > testerT[perm!=0]/(testerT[perm!=0]+testerF[perm!=0]))),
                  list(chr, inv, focalCluster, testCluster, mod, variable, thr)]

### save
  save(m.ag, m.ag.ag, file="/project/berglandlab/alan/environmental_ombibus_global/crossCluster_enrichment.Rdata")
