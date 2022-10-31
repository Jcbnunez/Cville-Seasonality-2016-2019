### libraries
#library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(5)
#library(ggplot2)

#goal- regather gwas and report number of sig snps and gif accross all chromsoomes. 
#### ---> These folders are the outcomes of the GWAS -- see code #3 in the github
pheno.dir2 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Dariaphenotypes"
pheno.dir3 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/Dariaphenotypes"
pheno.dir4 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/Cynthiaphenotypes"
pheno.dir5 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Cynthiaphenotypes"
#peaks = load("/scratch/bal7cg/Deficiency-Line-confirmation/PEAKS_for_ANALYSIS.Rdata")

inversions = readRDS("inv.dt")


#fiddilin- how to get out start stop points from inversions
#start = unlist(inversions[chrom == "2L"]$start)
#modify inversions table to have just one start/stop point for 3R
fixed3r = c("3R",11750567, 29031297,"total3r")

inversions = inversions %>% add_row(chrom = fixed3r[1], start = 11750567, stop = 29031297, invName = "fixed3")
#remove old rows
inversions = inversions[-c(3:5),]
#pheno.dir6 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/JanGWAS"

#pheno.files <- list.files(path=c(pheno.dir1, pheno.dir2, pheno.dir3, pheno.dir4), all.files=T, full.names=T, recursive=T)
pheno.files <- list.files(path=c(pheno.dir2, pheno.dir5), all.files=T, full.names=T, recursive=T)


#job <- 0 ==> sanity check
args = commandArgs(trailingOnly=TRUE)
job=as.numeric(args[1])-1



### thermal model
glm.fn <- paste("/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
#glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0", ".Rdata", sep="")
load(glm.fn)
#compare structure to Alan's new model files
# file = "/project/berglandlab/alan/dest_glm_final_nested_qb_alternatives_AIC/dest_glm_final_nested_qb_alternatives_AIC/job1.Rdata"
#  example = load(file)
### R5 -> R6 DGRP conversion table
liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
liftover <- fread(liftover.fn)
liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### iterate through phenotypes
gwas.o <- foreach(pheno.i=pheno.files, .errorhandling="remove")%dopar%{
  message(pheno.i)
  pheno.i <- pheno.files[1]
  
  ### load GWAS object
  pheno <- fread(pheno.i)
  
  ### convert GWAS output to R6
  setkey(pheno, SNP)
  setkey(liftover, SNP)
  
  pheno <- merge(pheno, liftover)
  setnames(pheno, c("dm6_chr", "dm6_pos"), c("chr", "pos"))
  
  ### merge & cleanup
  setkey(pheno, chr, pos)
  setkey(glm.out, chr, pos)
  m.all <- merge(glm.out, pheno)
  
  ### define thresholds
  #thrs <- sapply(c(1:9), function(x) x*10^c(-5:-1))
  #thrs <- expand.grid(thrs)[,1]
  thrs <- c(0.01, 0.05)
  
  ### run test for both the "aveTemp+year_factor" & "year_factor" GLMs
  gwas.o <- foreach(mod.i=c("aveTemp+year_factor", "year_factor"))%do%{
    # mod.i <-"aveTemp+year_factor"
    m <- m.all[!is.na(rnp.clean) & mod==mod.i]
    m[,glm.rnp:=rank(p.lrt)/(length(p.lrt)+1)]
    m[,gwas.rnp:=rank(PVAL)/(length(PVAL)+1)]
    
    if(mod.i=="aveTemp+year_factor") m[,glm.beta:=as.numeric(tstrsplit(b, ";")[[5]])]
    if(mod.i=="year_factor") m[,glm.beta:=NA]
    
    
    ### enrichment
    en <- foreach(thr=thrs, .combine="rbind")%do%{
      m[,
        list(TT=sum(glm.rnp<=thr & gwas.rnp<=thr),
             TF=sum(glm.rnp>=thr & gwas.rnp<=thr),
             FT=sum(glm.rnp<=thr & gwas.rnp>=thr),
             FF=sum(glm.rnp>=thr & gwas.rnp>=thr),
             thr=thr),
        list(inv=(invName!="none"), chr)]
    }
    
    
    en[,or:=(TT/TF)/(FT/FF)]
    en[,fet.p:=apply(en, 1, function(x) fisher.test(matrix(as.numeric(c(x[3], x[4], x[5], x[6])), byrow=T, nrow=2))$p.value)]
    # en[inv==T & chr=="2L"][order(thr)]
    
    
    ### sign test
    st <- foreach(thr=thrs, .combine="rbind")%do%{
      m[glm.rnp<=thr & gwas.rnp<=thr,
        list(st.T=sum(sign(glm.beta)==sign(SCORE), na.rm=T),
             st.F=sum(sign(glm.beta)!=sign(SCORE), na.rm=T),
             thr=thr),
        list(inv=(invName!="none"), chr)]
    }
    st[,prop:=(st.T)/(st.T+st.F)]
    # st[inv==T & chr=="2L"][order(thr)]
    
    ### single GLM or GWAS enrichments
    sen <- foreach(thr=thrs, .combine="rbind")%do%{
      m[,
        list(glm_T=sum(glm.rnp<=thr),
             glm_F=sum(glm.rnp>=thr),
             gwas_T=sum(gwas.rnp<=thr),
             gwas_F=sum(gwas.rnp>=thr),
             thr=thr),
        list(inv=(invName!="none"), chr)]
    }
    sen[,glm_prop:=glm_T/(glm_T+glm_F)]
    sen[,gwas_prop:=gwas_T/(gwas_T+glm_F)]
    
    
    
    ### merge output
    setkey(en, inv, chr, thr)
    setkey(st, inv, chr, thr)
    setkey(sen, inv, chr, thr)
    
    gwas.o <- merge(en, st)
    gwas.o <- merge(gwas.o, sen)
    
    gwas.o[,gwas.pheno := tstrsplit(pheno.i, "/") %>% last %>% gsub(".txt", "", .)]
    gwas.o[,glm.mod:=mod.i]
    
    ### return
    return(gwas.o)
  }
  gwas.o <- rbindlist(gwas.o)
  
  #gwas.o[inv==T & chr=="2L"][order(thr)]
  #gwas.o[,glm.pop:="VA_ch"]
  gwas.o[,glm.perm:=job]
  
  ### return
  return(gwas.o)
}
gwas.o <- rbindlist(gwas.o)

### save
system("mkdir gwas_glm_merge")
save(gwas.o, file=paste("./gwas_glm_merge/VA_ch_", job, ".Rdata", sep=""))
