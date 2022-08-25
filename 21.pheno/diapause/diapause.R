### libraries
  library(data.table)
  library(foreach)

### load Priscilla's data
  load("../gwas_top1percent.Rdat")

### R5 -> R6 DGRP conversion table
  liftover.fn <- "/project/berglandlab/Dmel_genomic_resources/liftOver_files/dest.all.PoolSNP.001.50.dm3.dm6.csv"
  liftover <- fread(liftover.fn)
  liftover[,SNP:=paste(dm3_chr, dm3_pos, "SNP", sep="_")]

### do liftover
  setnames(y, c("chr", "pos"), c("dm3_chr", "dm3_pos"))
  setkey(y, dm3_chr, dm3_pos)
  setkey(liftover, dm3_chr, dm3_pos)

  yl <- merge(y, liftover)
  setnames(yl, c("dm6_chr", "dm6_pos"), c("chr", "pos"))

### load thermal GLM
  job <- 0
  glm.fn <- paste("/project/berglandlab/thermal_glm_dest/dest_glm_final_nested_qb/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
  load(glm.fn)

  setkey(glm.out, chr, pos)

### merge
  setkey(yl, chr, pos)
  m <- merge(glm.out, yl, all.x=T)

### test
  m.ag <- m[chr=="2L"][mod=="aveTemp+year_factor"][pos>=16645762 & pos<=16955762][,
             list(TT=sum(!is.na(lod.geno) & rnp.clean<=0.05, na.rm=T),
                  TF=sum(!is.na(lod.geno) & rnp.clean>0.05,na.rm=T),
                  FT=sum(is.na(lod.geno) & rnp.clean<=0.05, na.rm=T),
                  FF=sum(is.na(lod.geno) & rnp.clean>0.05, na.rm=T)),
              list(mod, chr, inv=invName!="none", model, pop, phenotype)]

m.ag[perm.y==0][chr=="2L"][mod=="aveTemp+year_factor"][pop=="both"][phenotype=="diapause.bin9"]
