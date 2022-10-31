
# ijob -A berglandlab_standard -c10 -p standard --mem=40G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### get job info
jobId <- 15
args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])
#jobId = 1
### libraries
library(data.table)
library(foreach)
library(tidyr)
library(doMC)
registerDoMC(5)

### load in jobs files
jobs <- fread("/scratch/bal7cg/Deficiency-Line-confirmation/nogrm.newmodels.job.csv")
load( "/scratch/bal7cg/Deficiency-Line-confirmation/snp_dt_25percMissing.Rdata")
setkey(snp.dt, id)

use <- apply(snp.dt[,c("VA_ch"), with=F],
             1, any)
snp.dt <- snp.dt[use]
colnames(snp.dt)[1] = "variant.id"
setkey(snp.dt, variant.id)
snp.dt = snp.dt[,c(1:3)]
### load in  data
### phenotype files


pheno.files <- jobs[job==jobId]$gwas
# jobs[grepl("ActivityLevel_Standard-BasalActivity_F.nogrms.jan", gwas)]

#### thermal model
#glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
file = (jobs[job==jobId]$glm[1])
glm.out = readRDS(file)
#use snpid file to correctly merge in pos and chr
setkey(snp.dt, variant.id)
setkey(glm.out, variant.id)
glm.out = merge(glm.out,snp.dt)
### R5 -> R6 DGRP conversion table
liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
liftover <- fread(liftover.fn)
liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### iterate through phenotypes
gwas.win.o <- foreach(pheno.i=pheno.files, .errorhandling="remove")%do%{
  # pheno.i  = pheno.files[1]
  message(paste(which(pheno.i==pheno.files), " / ", length(pheno.files)))
  # pheno.i <- pheno.files[grepl("ActivityLevel_Standard-BasalActivity_F.nogrms.jan", pheno.files)]
  
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
  
  ### define windows
  
  ### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  setkey(m.all, "chr")
  wins <- foreach(chr.i=c("2L","2R","3L","3R"), .combine="rbind", .errorhandling="remove")%do%{
    #chr.i = "3R"
    tmp <- m.all[J(chr.i)]
    data.table(chr=chr.i,
               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)
  
  
  ### run test for both the "aveTemp+year_factor" & "year_factor" GLMs
  gwas.win.o <- foreach(mod.x=c(3,4))%do%{
    
    # mod.x = 3
    m <- m.all[ mod.i==mod.x]
    m[,glm.rnp:=rank(p_lrt)/(length(p_lrt)+1)]
    m[,gwas.rnp:=rank(PVAL)/(length(PVAL)+1)]
    
    
    setkey(m, chr, pos)
    
    gwas.win.o <- foreach(win.i=1:dim(wins)[1], .errorhandling="remove")%dopar%{
      if(win.i%%100==0) message(win.i)
      # win.i <- 103
      thr <- 0.05
      #filter out to data within the window
      tmp = m[chr == wins[win.i]$chr][pos >= wins[win.i]$start][pos <= wins[win.i]$end]
      
      ### enrichment
      en <- tmp[,
                list(TT=sum(glm.rnp<=thr & gwas.rnp<=thr),
                     TF=sum(glm.rnp>=thr & gwas.rnp<=thr),
                     FT=sum(glm.rnp<=thr & gwas.rnp>=thr),
                     FF=sum(glm.rnp>=thr & gwas.rnp>=thr),
                     thr=thr),
                list( chr)]
      en[,or:=(TT/TF)/(FT/FF)]
      en[,fet.p:=fisher.test(matrix(c(TT, TF, FT, FF), byrow=T, nrow=2))$p.value]
      
      # en[inv==T & chr=="2L"][order(thr)]
      
      
      ### sign test
      st <- tmp[glm.rnp<=thr & gwas.rnp<=thr,
                list(st.T=sum(sign(b_temp)==sign(SCORE), na.rm=T),
                     st.F=sum(sign(b_temp)!=sign(SCORE), na.rm=T),
                     thr=thr),
                list(chr)]
      st[,prop:=(st.T)/(st.T+st.F)]
      st[,prop.p:=binom.test(st.T, st.T+st.F, .5)$p.value]
      # st[inv==T & chr=="2L"][order(thr)]
      
      ### merge output
      setkey(en,  chr, thr)
      setkey(st,  chr, thr)
      
      gwas.o <- merge(en, st)
      
      ###Adam modification area###
      #create a vector list of TT snps
      tmp$TT.status = ifelse(tmp$glm.rnp <= thr & tmp$gwas.rnp <= thr, T, F)
      TTsnps = as.vector(tmp[TT.status == T]$pos)
      gwas.o$snp.positions = toString(TTsnps)
      #############################################
      
      gwas.o[,gwas.pheno := tstrsplit(pheno.i, "/") %>% last %>% gsub(".txt", "", .)]
      gwas.o[,glm.mod:=mod.x]
      gwas.o[,start:=wins[win.i]$start]
      gwas.o[,end:=wins[win.i]$end]
      gwas.o[,win.i:=win.i]
      
      ### return
      return(gwas.o)
    }
    gwas.win.o <- rbindlist(gwas.win.o)
    
    #gwas.o[inv==T & chr=="2L"][order(thr)]
    # gwas.win.o[,glm.pop:="VA_ch"]
    gwas.win.o[,glm.perm:=jobs[job==jobId]$glm[1]]
    
    ### return
    return(gwas.win.o)
  }
  gwas.win.o <- rbindlist(gwas.win.o)
  
}
gwas.win.o <- rbindlist(gwas.win.o)

warnings()

### save
save(gwas.win.o, file=paste("/scratch/bal7cg/Deficiency-Line-confirmation/newmodel.sliding.data/withoutgrm_window_", jobId, ".Rdata", sep=""))
