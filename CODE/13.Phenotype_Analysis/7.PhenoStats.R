### libraries
#library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(4)
#library(ggplot2)
#this script gathers together gwas files




### PATH is the path two the gmmat files to be gathered
PATH = "/scratch/bal7cg/Yang_Adam/GWAS_withoutGRMs/PermsOutputsWithoutGRMs"
#load in the taskid table used to create this set of gmmat files
task.id = fread("/scratch/bal7cg/Yang_Adam/perm.gwas.ref.table.txt")


### define significant p values to look at
p <- 0.00001
#run a set of nested foreach loops. the outer loop opens finds and opens each gmmat file. The inner loop runs from chromosome to chromosome, and compiles a set of summary statistics including the genomic inflation factor, the Pvalue and location of the most significant snp, and the # of snps passing threshold and # of snps total. These values are bound together in a list, with missing values removed
o.statistics <- foreach(i = c(1:1310),.combine = "rbind", .errorhandling = "remove")%dopar%{
  i = 2
  print(i)
  
  phenotypes = task.id[array.id == id]$pheno
  permutations = task.id[array.id == id]$perms
  
  ### change to whatever file naming format the files are named
  #file should match the "output" variable from gmmat generating script
  file <- paste0(PATH,"/", phenotypes, ".perms.nogrms/",phenotypes,".nogrms.perms.", permutations, ".txt")
  
  #load in score file
  gmmat = fread(file)
  gmmat <- gmmat[is.na(gmmat$PVAL) == F][,c("SNP","CHR","POS","AF","PVAL"),with=F] #remove missing data
  #clean previous assignmations to "gmmat" to conserve mememory
  gc(gmmat)
  print(paste(phenotypes, dim(gmmat)))
  
  
  chromout <- foreach(chr = unique(c("2L","2R","3L","3R","X")),.combine = "rbind") %do% {
    
    #filter down to a certain chromsomome
    gmmat.chr = gmmat[CHR == chr]
    #find genomic inflation
    chisq <- qchisq(1-gmmat.chr$PVAL,1)
    GE = median(chisq)/qchisq(0.5,1)
    #find top snp
    #gmmat.chr = gmmat.chr %>% arrange(PVAL)
    minpval = min(gmmat.chr$PVAL)
    top.id = gmmat.chr[PVAL == minpval]$SNP
    
    #find # of snps with pval < 0.001
    sig.number = dim(gmmat.chr[PVAL < p]) [1]
    #find # of snps total
    totalsnpnumber = dim(gmmat.chr)[1]
    #put together summary statistics
    values = data.table (
      chromosome = chr ,
      #phenotype = phenotype.tag,
      GIF = GE,
      Top.id = top.id,
      Top.pval = minpval,
      Sig.number = sig.number,
      Sig.threshold = p,
      Totalsnpnumber = totalsnpnumber
    )
    values
  }
  chromout[,phenotype:=phenotypes]
  chromout[,permutation:=permutations]
  
  chromout
  
}
#save summarized data into a table
write.table(o.statistics,file=paste0("/scratch/bal7cg/Yang_Adam/gathers/novperm.stats.withoutgrm.txt"),quote=F,row.names=F,col.names=T,sep="\t")



