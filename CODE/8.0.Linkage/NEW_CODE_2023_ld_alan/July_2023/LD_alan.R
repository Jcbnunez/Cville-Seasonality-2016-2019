# ijob -A berglandlab_standard -c20 -p standard --mem=30G
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(SNPRelate)
  library(data.table)

### convert to snprelate GDS
  vcf.fn <- "/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz"
  gds.fn <- "/scratch/aob2x/M_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.snprelate.gds"

  snpgdsVCF2GDS(vcf.fn, gds.fn, verbose=T)

### open
  genofile <- snpgdsOpen(gds.fn)

### get SNP list

  snp.dt <- as.data.table(snpgdsSNPList(genofile, sample.id=NULL))
  snp.dt[chromosome=="2L"][afreq>.05 & afreq<.95]

### make sliding window LD matrix
ldmat <- snpgdsLDMat(genofile, 
            sample.id=NULL,
            snp.id=snp.dt[chromosome=="2L"][afreq>.05 & afreq<.95]$snp.id[1:1000],
            slide=50,
            method=c("composite"), mat.trim=T,
            num.thread=20L, with.id=TRUE, verbose=TRUE)
str(ldmat)



###
scp aob2x@rivanna.hpc.virginia.edu:~/ldmat.Rdata ~/.

library(ggplot2)

load("~/ldmat.Rdata")
image(t(ldmat$LD), col=terrain.colors(16))
  snpgdsLDpair(snp1=)
