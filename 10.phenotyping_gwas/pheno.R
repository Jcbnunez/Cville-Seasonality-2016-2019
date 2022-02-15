# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### get job info
  # job <- 1
  args = commandArgs(trailingOnly=TRUE)
  job=as.numeric(args[1])-1

### libraries
  library(data.table)
  library(foreach)
  library(tidyr)
  library(doMC)
  registerDoMC(10)

### load in  data
  ### phenotype files
    #pheno.fn <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS/HighThermalToleranceExtreme_VaryingWithTemperature_F.original.nogrms/HighThermalToleranceExtreme_VaryingWithTemperature_F.nogrms.original.txt"
    pheno.dir1 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS"
    pheno.dir2 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/OriginalOutputsWithoutGRMS"
    pheno.dir3 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/NewPhenoOutputsWithGRMS"
    pheno.dir4 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/OriginalOutputsWithGRMS"

    pheno.files <- list.files(path=c(pheno.dir1, pheno.dir2, pheno.dir3, pheno.dir4), all.files=T, full.names=T, recursive=T)

  ### thermal model
    glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
    load(glm.fn)

  ### R5 -> R6 DGRP conversion table
    liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
    liftover <- fread(liftover.fn)
    liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### iterate through phenotypes
  gwas.o <- foreach(pheno.i=pheno.files, .errorhandling="remove")%dopar%{
    message(pheno.i)
    # pheno.i <- pheno.files[74]

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
      thrs <- sapply(c(1:9), function(x) x*10^c(-5:-1))
      thrs <- expand.grid(thrs)[,1]


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

        ### merge output
          setkey(en, inv, chr, thr)
          setkey(st, inv, chr, thr)

          gwas.o <- merge(en, st)

          gwas.o[,gwas.pheno := tstrsplit(pheno.i, "/") %>% last %>% gsub(".txt", "", .)]
          gwas.o[,glm.mod:=mod.i]

        ### return
          return(gwas.o)
      }
      gwas.o <- rbindlist(gwas.o)

      #gwas.o[inv==T & chr=="2L"][order(thr)]
      gwas.o[,glm.pop:="VA_ch"]
      gwas.o[,glm.perm:=job]

    ### return
      return(gwas.o)
  }
  gwas.o <- rbindlist(gwas.o)

### save
  save(gwas.o, file=paste("/scratch/aob2x/gwas_glm_merge/VA_ch_", job, ".Rdata", sep=""))
