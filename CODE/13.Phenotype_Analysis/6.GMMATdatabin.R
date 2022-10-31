
#this script takes the gmmmat gather files and gathers them again into one data file
#with addition columns for observed/permutation and GRM status
#set wd to location of gather files
setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/November_objects/")
library(tidyverse)
library(data.table)
ob.nogrms = fread("./novperm.stats.withoutgrm.txt")
ob.grms = fread("./novperm.stats.withgrm.txt")

#ob.nogrms = ob.nogrms[,c("phenotype","chromosome","Top.pval","Sig.number","GIF", "permutation"),with=F]
ob.nogrms[,GRMs:=F]
#setnames(ob.grms,"GE","GIF")
#ob.grms = ob.grms[,c("phenotype","chromosome","Top.pval","Sig.number","GIF"),with=F]
ob.grms[,GRMs:=T]
setkey(ob.nogrms,phenotype,chromosome, permutation)
setkey(ob.grms,phenotype,chromosome, permutation)
perm.data = merge(ob.nogrms,ob.grms)

perm.data[,delta.GIF := GIF.x-GIF.y]
perm.data[,delta.sig := (Sig.number.x / Totalsnpnumber.x) - (Sig.number.y / Totalsnpnumber.y)]
#change permutation column to a simple "permutation"
perm.data$permutation = as.character(perm.data$permutation)
perm.data$permutation = "permutation"
#lets add in the old and new pheno (observed)
old.pheno.observed.nogrm = fread("pheno.stats.withoutgrm.txt")
old.pheno.observed.grm = fread("pheno.stats.withgrm.txt")
new.pheno.observed.nogrm = fread("newpheno.stats.withoutgrm.txt")
new.pheno.observed.grm = fread("pheno.stats.withgrm.txt")
#create variable for GRM status
old.pheno.observed.grm$GRMs = T
old.pheno.observed.nogrm$GRMs = F
new.pheno.observed.grm$GRMs = T
new.pheno.observed.nogrm$GRMs = F
#rbind data for grm and no grm together
observed.grm = rbind(old.pheno.observed.grm, new.pheno.observed.grm)
observed.nogrm = rbind(old.pheno.observed.nogrm, new.pheno.observed.nogrm)
#merge in same way as permutated data
setkey(observed.grm,phenotype,chromosome)
setkey(observed.nogrm,phenotype,chromosome)
observed.data = merge(observed.nogrm, observed.grm)
#create columns for delta gif and delta sig and permutation (false)
observed.data[,delta.GIF := GIF.x-GIF.y]
observed.data[,delta.sig := (Sig.number.x / Totalsnpnumber.x) - (Sig.number.y / Totalsnpnumber.y)]
observed.data[,permutation := "observed"]
#filer out phenotypes we're not intersted in
updatednames = readRDS("updatedphenonames")
updatednames = data.table(
  phenotype = updatednames
  
)
observed.data = merge(updatednames, observed.data, by = "phenotype")
#finally, rbind data together
all.data = rbind(perm.data, observed.data)
saveRDS(all.data, "allgatheredgwasstats")
ob.order = ob.all[order(ob.all[chromosome == "2L"]$delta.GIF)]
