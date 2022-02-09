### Run the Berry and Kreirman test
### module load intel/18.0 intelmpi/18.0 R/3.6.3; R
args = commandArgs(trailingOnly=TRUE)
jobId=as.numeric(args[1])


## jobId=1

### libraries
library(data.table)
library(gdata)
library(lubridate)
library(foreach)
library(SeqArray)
library(doMC)
registerDoMC(2)
library(tidyr)

### load in Inversion Marker Sets
inv.mark <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)


inv.mark[,chr:=tstrsplit(V1, "_")[[1]]]
inv.mark[,pos:=as.numeric(tstrsplit(V1, "_")[[2]])]

setkey(inv.mark, chr, pos)
inv.mark[,inv:="In2Lt"]

### load in temperature GLM model

root_path <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out."
### Extratc SNPs
iterations = c(0,1:98,100)
print(iterations[i])

pop= "VA_ch"
i=1
model = "aveTemp+year_factor"
p_tresh=0.05

load(paste(root_path,
           pop,
           "_", 
           iterations[i], 
           ".Rdata", sep = ""))

category = ifelse(iterations[i] == 0, "Obs", "Per")
print(category)

#load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")
setkey(glm.out, chr, pos)

glm.out <- merge(glm.out, inv.mark, all.x=T)
glm.out[is.na(inv), inv:="none"]
glm.out[,V1:=paste(chr, pos, "SNP", sep="_")]

table(glm.out$rnp.clean<=p_tresh, glm.out$inv)

glm.out.short <- glm.out[mod==model][chr=="2L"][rnp.clean<=p_tresh | inv!="none"]

table(glm.out$inv)
table(glm.out.short$inv)

### open PoolSeq GDS
pool.gds <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### PoolSeq SNP.dt
load("./snp_dt_25percMissing.Rdata")
pool.snp.dt <- snp.dt
#rm(snp.dt)
#pool.snp.dt[,V1:=paste(chr, pos, "SNP", sep="_")]

### open Phased Individual Seq GDS
ind.gds <- seqOpen("/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.gds", allow.duplicate=T)

### individual seq SNP tables
ind.snp.dt <- data.table(chr=seqGetData(ind.gds, "chromosome"),
                         pos=seqGetData(ind.gds, "position"),
                         variant.id=seqGetData(ind.gds, "variant.id"))
#ind.snp.dt[,V1:=paste(chr, pos, "SNP", "_")]

setkey(ind.snp.dt,  chr, pos)
setkey(inv.mark,    chr, pos)
setkey(pool.snp.dt, chr, pos)
setkey(glm.out.short, chr, pos)

glm.out.short <- merge(glm.out.short, ind.snp.dt)


### weather data
load("./weatherAve.Rdata")
setnames(weather.ave, "V1", "sampleId")

### samps
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

samps[locality=="UA_od", locality:="UA_Ode"]

samps <- merge(samps, weather.ave[,c("aveTemp", "sampleId")], by="sampleId")

### subset GLM.out to SNPs in individual genotyping
#  setkey(glm.out, chr, pos)
#  setkey(ind.snp.dt, chr, pos)
#
#  merge(glm.out, ind.snp.dt[,c("chr", "pos", "variant.id"), with=F], no.match="NA")


### function to get observed allele frequencies
getAF <- function(sampleId, variantId) {
  
  seqSetFilter(pool.gds,
               sample.id=sampleId,
               variant.id=variantId)
  
  ad <- seqGetData(pool.gds, "annotation/format/AD")
  dp <- seqGetData(pool.gds, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp)[,1],
                   sampleId=rep(seqGetData(pool.gds, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(pool.gds, "variant.id"), each=dim(ad$data)[1]))
  
  afi <- merge(af, pool.snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  return(afis[,c("sampleId", "variant.id", "chr", "pos", "aveTemp", "nEff", "af_nEff"), with=F])
  
}

### function to get expected frequecies
getExpectedFrequencies <- function(driver, passenger, af.driver=NULL) {
  ### driver <- "2L_13185689_SNP"; passenger <- "2L_13186585_SNP"
  ### driver <- "2L_13185689_SNP"; passenger <- "2L_13185689_SNP"
  setkey(driver, chr, pos)
  setkey(passenger, chr, pos)
  
  ### subset filter
  seqSetFilter(ind.gds,
               variant.id=ind.snp.dt[J(driver)]$variant.id)
  locus1 <- seqGetData(ind.gds, "genotype")[,,1] %>%  t()
  
  seqSetFilter(ind.gds,
               variant.id=ind.snp.dt[J(passenger)]$variant.id)
  locus2 <- seqGetData(ind.gds, "genotype")[,,1] %>%  t()
  
  ### generate LD table
  tab.s1 <- table(driver=locus1[,1]==1, passenger=locus2[,1]==1)
  tab.s2 <- table(driver=locus1[,2]==1, passenger=locus2[,2]==1)
  
  tab <- tab.s1 + tab.s2
  
  ldtab <- tab/sum(tab)  ### FALSE==REF; TRUE==ALT
  p1=ldtab[1,1] + ldtab[1,2]
  p2=ldtab[2,1] + ldtab[2,2]
  q1=ldtab[1,1] + ldtab[2,1]
  q2=ldtab[1,2] + ldtab[2,2]
  
  D=ldtab[1,1] -p1*q1; D
  r2=D^2/(p1*q1*p2*q2); r2
  
  ### get observed allele frequencies
  message("Driver Allele Freqs")
  
  if(is.null(af.driver)) {
    af.driver <- getAF(sampleId=samps[locality=="VA_ch"][year>=2014]$sampleId,
                       variantId=pool.snp.dt[J(driver)]$id)
  }
  
  af.passenger <- getAF(sampleId=samps[locality=="VA_ch"][year>=2014]$sampleId,
                        variantId=pool.snp.dt[J(passenger)]$id)
  
  ### get expected passenger allele frequencies
  af.passenger[,exp.af:= (ldtab[2,2]/sum(ldtab[,2]))*(af.driver$af_nEff) +
                 (ldtab[1,2]/sum(ldtab[,2]))*(1-af.driver$af_nEff)]
  
  af.passenger[,exp.af_nEff:=round(exp.af*nEff)/nEff]
  
  af.passenger[,driver:=driver$V1]
  af.passenger[,passenger:=passenger$V1]
  
  ### tack in summary stat
  af.passenger[,D:=D]
  af.passenger[,r2:=r2]
  
  ### return
  return(af.passenger)
}

### function to run GLM on binomial resampling of expected (or observed) allele frequencies
simulateGLM <- function(dat, nBoot=100) {
  #
  
  dat <- dat[!is.na(exp.af) & !is.na(nEff)]
  t1.obs <- glm(af_nEff~aveTemp, data=dat, family=binomial(), weights=dat$nEff)
  #t1.obs <- lm(asin(sqrt(af_nEff))~aveTemp, data=dat)
  
  
  sim.out <- foreach(i=1:nBoot, .combine="rbind")%dopar%{
    dat[,sim.af_nEff:=rbinom(length(dat$nEff), dat$nEff, dat$exp.af)/dat$nEff]
    
    t1.sim <- glm(sim.af_nEff~aveTemp, data=dat, family=binomial(), weights=dat$nEff)
    #t1.sim <- lm(asin(sqrt(sim.af_nEff))~aveTemp, data=dat)
    
    data.table(boot=i, beta=coef(t1.sim)[2], pseudoR2=(1 - (t1.sim$deviance/t1.sim$null.deviance)))
    #data.table(boot=i, beta=coef(t1.sim)[2], pseudoR2=summary(t1.sim)$r.squared)
  }
  
  sim.out <- rbind(sim.out, data.table(boot=0, beta=coef(t1.obs)[2], pseudoR2=(1 - (t1.obs$deviance/t1.obs$null.deviance))))
  #sim.out <- rbind(sim.out, data.table(boot=0, beta=coef(t1.obs)[2], pseudoR2=summary(t1.obs)$r.squared))
  
  return(sim.out)
}

### LM simulation
simulateLM <- function(dat, nBoot=100) {
  #
  
  dat <- dat[!is.na(exp.af) & !is.na(nEff)]
  #t1.obs <- glm(af_nEff~aveTemp, data=dat, family=binomial(), weights=dat$nEff)
  t1.obs <- lm(asin(sqrt(af_nEff))~aveTemp, data=dat)
  
  
  sim.out <- foreach(i=1:nBoot, .combine="rbind")%dopar%{
    dat[,sim.af_nEff:=rbinom(length(dat$nEff), dat$nEff, dat$exp.af)/dat$nEff]
    
    #t1.sim <- glm(sim.af_nEff~aveTemp, data=dat, family=binomial(), weights=dat$nEff)
    t1.sim <- lm(asin(sqrt(sim.af_nEff))~aveTemp, data=dat)
    
    #data.table(boot=i, beta=coef(t1.sim)[2], pseudoR2=(1 - (t1.sim$deviance/t1.sim$null.deviance)), class="sim")
    data.table(boot=i, beta=coef(t1.sim)[2], pseudoR2=summary(t1.sim)$r.squared)
  }
  
  #sim.out <- rbind(sim.out, data.table(boot=0, beta=coef(t1.obs)[2], pseudoR2=(1 - (t1.obs$deviance/t1.obs$null.deviance))))
  sim.out <- rbind(sim.out, data.table(boot=0, beta=coef(t1.obs)[2], pseudoR2=summary(t1.obs)$r.squared))
  
  return(sim.out)
}

### what about just a correlation of expected and observed frequencies at the passenger site?
exp_obs_cor <- function(dat, nBoot=100) {
  dat <- dat[!is.na(exp.af) & !is.na(nEff)]
  ct <- cor.test(asin(sqrt(dat$exp.af_nEff)), asin(sqrt(dat$af_nEff)))
  data.table(cor=ct$estimate, p=ct$p.value)
}

### likelihood function
ll <- function(dat) {
  dat.tmp <- dat[exp.af_nEff!=0 & exp.af_nEff!=1]
  data.table(ll_exp=sum(dbinom(dat.tmp$af_nEff*dat.tmp$nEff, dat.tmp$nEff, dat.tmp$exp.af_nEff, log=T), na.rm=T),
             ll_obs=sum(dbinom(dat.tmp$af_nEff*dat.tmp$nEff, dat.tmp$nEff, dat.tmp$af_nEff, log=T), na.rm=T),
             N=dim(dat.tmp)[1])
  
}

### define driver and passenger sets
drivers <- rbind(
  data.table(chr=glm.out.short[r.clean<=250]$chr,
             pos=glm.out.short[r.clean<=250]$pos,
             V1=glm.out.short[r.clean<=250]$V1,
             rnp.clean=glm.out.short[r.clean<=250]$rnp.clean,
             rank.clean=glm.out.short[r.clean<=250]$r.clean,
             class="aveTemp", key="chr,pos"),
  
  data.table(chr=glm.out.short[inv!="none"]$chr,
             pos=glm.out.short[inv!="none"]$pos,
             V1=glm.out.short[inv!="none"]$V1,
             rnp.clean=glm.out.short[inv!="none"]$rnp.clean,
             rank.clean=glm.out.short[inv!="none"]$r.clean,
             class="inv", key="chr,pos"))

dim(drivers)

passengers <- data.table(chr=glm.out.short$chr,
                         pos=glm.out.short$pos,
                         V1=glm.out.short$V1,
                         rnp.clean= glm.out.short$rnp.clean,
                         rank.clean=glm.out.short$r.clean,
                         class="bg", key="chr,pos")

### iterate through
setkey(passengers, chr, pos)



#### RUN THE DOPAR THE BK TEST
#### 
bk.o <- foreach(driver.i=c(1:dim(drivers)[1])[jobId], .combine="rbind", .errorhandling="remove")%do%{
  #driver.i=1; passenger.i<-12
  
  af.driver <- getAF(sampleId=samps[locality=="VA_ch"][year>=2014]$sampleId,
                     variantId = pool.snp.dt[J(drivers[driver.i])]$id)
  
  foreach(passenger.i=c(1:dim(passengers)[1]), .combine="rbind", .errorhandling="remove")%dopar%{
    
    message(paste(driver.i, passenger.i, sep=" / "))
    
    dat.pd <- getExpectedFrequencies(driver=drivers[driver.i], passenger=passengers[passenger.i], af.driver=af.driver)
    sim.glm <- simulateGLM(dat=dat.pd, nBoot=100) ## GLM version
    sim.glm[,mod:="glm"]
    
    sim.lm <- simulateLM(dat=dat.pd, nBoot=100) ## GLM version
    sim.lm[,mod:="lm"]
    sim <- rbind(sim.glm, sim.lm)
    
    sim.ag <- sim[,list(beta.pr=mean(beta[boot==0]>beta[boot!=0]),
                        beta.obs=beta[boot==0],
                        beta.sim.mu=mean(beta[boot!=0]), beta.sim.sd=sd(beta[boot!=0]),
                        pseudoR2.pr=mean(pseudoR2[boot==0]>pseudoR2[boot!=0]),
                        pseudoR2.obs=pseudoR2[boot==0],
                        pseudoR2.sim.mu=mean(pseudoR2[boot!=0]), pseudoR2.sim.sd=sd(pseudoR2[boot!=0]),
                        passenger=passengers[passenger.i]$V1, driver=drivers[driver.i]$V1,
                        passenger.class=passengers[passenger.i]$class, driver.class=drivers[driver.i]$class),
                  list(mod)]
    
    llo <- ll(dat=dat.pd)
    sim.ag <- cbind(sim.ag, llo)
    sim.ag[,D:=dat.pd$D[1]]
    sim.ag[,r2.ld:=dat.pd$r2[1]]
    sim.ag[,pai:=passenger.i]
    sim.ag[,di:=driver.i]
    
    sim.ag[,driver.rnp:=drivers[driver.i]$rnp.clean]
    sim.ag[,driver.rank:=drivers[driver.i]$rank.clean]
    
    sim.ag[,passenger.rnp:= passengers[passenger.i]$rnp.clean]
    sim.ag[,passenger.rank:=passengers[passenger.i]$rank.clean]
    
    return(sim.ag)
  }
  
} ### close the do loop

head(bk.o)
bk.o[,chisq:= 2*(ll_exp/ll_obs)]
bk.o[,p:=1 - pchisq(chisq, 2)]

### save
message("saving file")
save(bk.o, file=paste("/scratch/aob2x/bk_out/bk_", jobId, ".Rdata", sep=""))
