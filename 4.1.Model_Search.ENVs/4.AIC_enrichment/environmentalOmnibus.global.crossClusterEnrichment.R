# ijob -A berglandlab_standard -c20 -p standard --mem=50G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])


### libraries
  library(data.table)
  #library(gdata)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### get files
  fl <-  list.files("/project/berglandlab/alan/environmental_ombibus_global", full.names=T, recursive=T)
  fl <- as.data.table(fl[grepl("glmRNP", fl)])
  fl[,mod:=as.numeric(tstrsplit(last(tstrsplit(fl$V1, "/")), ";")[[2]])]
  fl[,variable:=tstrsplit(last(tstrsplit(fl$V1, "/")), ";")[[1]]]
  fl[,cluster:=gsub(".glmRNP.Rdata", "", tstrsplit(last(tstrsplit(fl$V1, "/")), ";")[[3]])]
  setkey(fl, mod, variable, cluster)

### which ones aren't there?
  load("/project/berglandlab/alan/environmental_ombibus_global/mod_var.Rdata")
  setkey(mod_var, mod, variable, cluster)
  mod_var[,job:=1:dim(mod_var)[1]]

  fl <- merge(fl, mod_var, all=T)
  table(fl[is.na(V1)]$mod, fl[is.na(V1)]$variable)
  fl[variable=="humidity.var" & mod==11 & is.na(V1)]

  fl <- fl[!is.na(V1)]

### get file paths
  setkey(fl, mod, variable)
  cville.fl <- fl[cluster=="5.Cville"][jobId]
  non_cville.fl <- fl[J(cville.fl)][cluster!="5.Cville"]

### load Cville
  load(file=cville.fl$V1)
  m1 <- glm.out

### iterate through
  m.ag <- foreach(i=1:dim(non_cville.fl)[1], .combine="rbind", .errorhandling="remove")%do%{
    fl.i <- non_cville.fl[i]
    message(fl.i$V1)
    load(fl.i$V1)
    m2 <- glm.out

    setkey(m1, chr, pos, perm, inv)
    setkey(m2, chr, pos, perm, inv)

    m <- merge(m1[,c("chr", "pos", "inv", "perm", "b_temp", "se_temp", "p_lrt", "rnp"), with=F],
               m2[,c("chr", "pos", "inv", "perm", "b_temp", "se_temp", "p_lrt", "rnp"), with=F])

    rm(m2, glm.out)

    ### summarize

      ### enrichment
        thrs <- c(0.01, 0.05)

        m.ag <- foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%do%{
          #thr.i <- thrs[1]

          m.ag <- m[,list(TT=sum(rnp.x<=thr.i & rnp.y<=thr.i),
                          TF=sum(rnp.x<=thr.i & rnp.y>thr.i),
                          FT=sum(rnp.x>thr.i & rnp.y<=thr.i),
                          FF=sum(rnp.x>thr.i & rnp.y>thr.i),
                          st.T=sum(sign(b_temp.x[rnp.x<=thr.i & rnp.y<=thr.i])==sign(b_temp.y[rnp.x<=thr.i  & rnp.y<=thr.i])),
                          st.F=sum(sign(b_temp.x[rnp.x<=thr.i & rnp.y<=thr.i])!=sign(b_temp.y[rnp.x<=thr.i  & rnp.y<=thr.i])),
                          focalT=sum(rnp.x<=thr.i),
                          focalF=sum(rnp.x>thr.i),
                          testerT=sum(rnp.y<=thr.i),
                          testerF=sum(rnp.x>thr.i),
                          thr=thr.i),
                      list(chr, inv, perm)]
          m.ag[,or:=(TT/TF)/(FT/FF)]
          m.ag[,st.pr:=st.T/(st.T+st.F)]
          m.ag[,focalCluster:=cville.fl$cluster]
          m.ag[,testCluster:=fl.i$cluster[]]
          m.ag[,mod:=cville.fl$mod]
          m.ag[,variable:=cville.fl$variable]
      }
      return(m.ag)
    }

### save
  mainDir <- "/scratch/aob2x/environmental_ombibus_global"
  outDir <- paste(mainDir, "crossClusterEnrichment", sep="/")
  dir.create(file.path(outDir), showWarnings = FALSE)

  save(m.ag, file=paste(outDir, "/job", jobId, ".Rdata", sep=""))

  message(paste(outDir, "/job", jobId, ".Rdata", sep=""))
