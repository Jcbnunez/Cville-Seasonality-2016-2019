# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
 library(data.table)
 library(gdata)
 library(lubridate)
 library(foreach)
 #library(SeqArray)
 #library(glmnet)
 library(doMC)
 registerDoMC(1)

### load in SNP file
  load("~/Overwintering_18_19/SNP_filtering/snp_dt.Rdata")

### define sample sets
  targetLocales <- list("FI_Aka"=list(locality="FI_Aka",   minYear=2000),
                      "DE_Bro"=list(locality="DE_Bro",   minYear=2000),
                      "VA_ch"=list(locality="VA_ch",     minYear=2014),
                      "PA_li"=list(locality="PA_li",     minYear=2000),
                      "DE_Mun"=list(locality="DE_Mun",   minYear=2000),
                      "UA_Ode"=list(locality=c("UA_Ode"), minYear=2010),
                      "TR_Yes"=list(locality="TR_Yes",   minYear=2000))

### job parameters
  nJobs <- 999
  nPerm <- 10

  
