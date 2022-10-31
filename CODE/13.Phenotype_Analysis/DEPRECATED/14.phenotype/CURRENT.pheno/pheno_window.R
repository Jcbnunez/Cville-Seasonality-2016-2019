# ijob -A berglandlab_standard -c10 -p standard --mem=40G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### get job info
  # jobId <- 15
  args = commandArgs(trailingOnly=TRUE)
  jobId=as.numeric(args[1])

### libraries
  library(data.table)
  library(foreach)
  library(tidyr)
  library(doMC)
  registerDoMC(1)

### load in jobs files
  jobs <- fread("~/Overwintering_18_19/pheno/pheno_jobs2.csv")

### load in  data
  ### phenotype files
    #pheno.fn <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS/HighThermalToleranceExtreme_VaryingWithTemperature_F.original.nogrms/HighThermalToleranceExtreme_VaryingWithTemperature_F.nogrms.original.txt"
    #pheno.dir1 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS"
    #pheno.dir2 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/OriginalOutputsWithoutGRMS"
    #pheno.dir3 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/NewPhenoOutputsWithGRMS"
    #pheno.dir4 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/OriginalOutputsWithGRMS"
    #
    #pheno.files <- list.files(path=c(pheno.dir1, pheno.dir2, pheno.dir3, pheno.dir4), all.files=T, full.names=T, recursive=T)
    #

    pheno.files <- jobs[job==jobId]$gwas
    # jobs[grepl("ActivityLevel_Standard-BasalActivity_F.nogrms.jan", gwas)]

  #### thermal model
    #glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
    load(jobs[job==jobId]$glm[1])

  ### R5 -> R6 DGRP conversion table
    liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
    liftover <- fread(liftover.fn)
    liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### iterate through phenotypes
  gwas.win.o <- foreach(pheno.i=pheno.files, .errorhandling="remove")%do%{
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
      wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%dopar%{
          tmp <- m.all[J(chr.i)]
          data.table(chr=chr.i,
                      start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                      end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
      }
      wins[,i:=1:dim(wins)[1]]
      dim(wins)


    ### run test for both the "aveTemp+year_factor" & "year_factor" GLMs
      gwas.win.o <- foreach(mod.i=c("aveTemp+year_factor"))%do%{

        # mod.i <-"aveTemp+year_factor"
        m <- m.all[!is.na(rnp.clean) & mod==mod.i]
        m[,glm.rnp:=rank(p.lrt)/(length(p.lrt)+1)]
        m[,gwas.rnp:=rank(PVAL)/(length(PVAL)+1)]

        if(mod.i=="aveTemp+year_factor") m[,glm.beta:=as.numeric(tstrsplit(b, ";")[[5]])]
        if(mod.i=="year_factor") m[,glm.beta:=NA]

        setkey(m, chr, pos)

        gwas.win.o <- foreach(win.i=1:dim(wins)[1], .errorhandling="remove")%dopar%{
          if(win.i%%100==0) message(win.i)
          # win.i <- 103
          thr <- 0.05

          tmp <- m[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

          ### enrichment
            en <- tmp[,
                  list(TT=sum(glm.rnp<=thr & gwas.rnp<=thr),
                     TF=sum(glm.rnp>=thr & gwas.rnp<=thr),
                     FT=sum(glm.rnp<=thr & gwas.rnp>=thr),
                     FF=sum(glm.rnp>=thr & gwas.rnp>=thr),
                     thr=thr),
                  list(inv=(invName!="none"), chr)]
            en[,or:=(TT/TF)/(FT/FF)]
            en[,fet.p:=fisher.test(matrix(c(TT, TF, FT, FF), byrow=T, nrow=2))$p.value]
            en <- en[!is.na(inv)]
            # en[inv==T & chr=="2L"][order(thr)]


          ### sign test
            st <- tmp[glm.rnp<=thr & gwas.rnp<=thr,
                    list(st.T=sum(sign(glm.beta)==sign(SCORE), na.rm=T),
                        st.F=sum(sign(glm.beta)!=sign(SCORE), na.rm=T),
                        thr=thr),
                    list(inv=(invName!="none"), chr)]
            st[,prop:=(st.T)/(st.T+st.F)]
            st[,prop.p:=binom.test(st.T, st.T+st.F, .5)$p.value]
            # st[inv==T & chr=="2L"][order(thr)]

          ### merge output
            setkey(en, inv, chr, thr)
            setkey(st, inv, chr, thr)

            gwas.o <- merge(en, st)

            gwas.o[,gwas.pheno := tstrsplit(pheno.i, "/") %>% last %>% gsub(".txt", "", .)]
            gwas.o[,glm.mod:=mod.i]
            gwas.o[,start:=wins[win.i]$start]
            gwas.o[,end:=wins[win.i]$end]
            gwas.o[,win.i:=win.i]

          ### return
            return(gwas.o)
        }
        gwas.win.o <- rbindlist(gwas.win.o)

        #gwas.o[inv==T & chr=="2L"][order(thr)]
        gwas.win.o[,glm.pop:="VA_ch"]
        gwas.win.o[,glm.perm:=jobs[job==jobId]$glm[1]]

        ### return
          return(gwas.win.o)
      }
      gwas.win.o <- rbindlist(gwas.win.o)

    }
  gwas.win.o <- rbindlist(gwas.win.o)

  warnings()

### save
  save(gwas.win.o, file=paste("/scratch/aob2x/gwas_glm_merge_window/run2_window_", jobId, ".Rdata", sep=""))
