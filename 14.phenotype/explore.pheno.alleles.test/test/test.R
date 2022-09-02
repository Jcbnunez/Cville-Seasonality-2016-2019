# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### get job info
  # job <- 0
  args = commandArgs(trailingOnly=TRUE)
  job=as.numeric(args[1])-1

### libraries
  library(data.table)
  library(foreach)
  library(tidyr)
  library(doMC)
  registerDoMC(5)

### load in  data
  ### phenotype files
    #pheno.fn <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS/HighThermalToleranceExtreme_VaryingWithTemperature_F.original.nogrms/HighThermalToleranceExtreme_VaryingWithTemperature_F.nogrms.original.txt"
    pheno.dir1 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS"
    pheno.dir2 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/OriginalOutputsWithoutGRMS"
    pheno.dir3 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/NewPhenoOutputsWithGRMS"
    pheno.dir4 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/OriginalOutputsWithGRMS"
    pheno.dir5 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/JanGWAS"
    pheno.dir6 <- "/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withGRMs/JanGWAS"

    pheno.files <- list.files(path=c(pheno.dir1, pheno.dir2, pheno.dir3, pheno.dir4, pheno.dir5, pheno.dir6), all.files=T, full.names=T, recursive=T)

  ### thermal model
    glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
    load(glm.fn)
