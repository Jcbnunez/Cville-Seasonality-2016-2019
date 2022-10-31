### libraries
library(data.table)
library(foreach)
library(tidyr)
#set working directory to directory for project
#setwd("/scratch/bal7cg/Yang_Adam")

#this script:1 creates directories to hold gwas files, 2: creates a reference table read by gwas scripts, 3: creates permuted phenotype data files for use in permutation gwas creation. 


### read in data files

INPUT= "."
#INPUT = "/scratch/bal7cg/Yang_Adam/"# path to working directory
#folder for observed GWAS with GRMs
ORIGINALOUTPUTGRMS = "./GWAS_withGRMs/NewPhenoOutputsWithGRMS/"
#folder for permutated GWAS with GRMs
PERMSOUTPUTGRMS = "./GWAS_withGRMs/PermsOutputsWithGRMS/"
#folder for observed GWAS without GRMs
ORIGINALOUTPUTNOGRMS = "./GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS/"
#folder for permutated GWAS without GRMs
PERMSOUTPUTNOGRMS = "./GWAS_withoutGRMs/PermsOutputsWithoutGRMs/"
#load in phenotype data
pheno.original = readRDS(paste0(INPUT,"updatedphenos"))
#fix names
#replace any / with _
#phenotype names need to be devoid of special characteristics and spaces
#these could cause directories to be created incorrectly
columns = colnames(pheno.original)
columns = gsub("/", "_", columns)
columns = gsub("\\(", "-", columns)
columns = gsub("\\)", "-", columns)
colnames(pheno.original) = columns

#example of fixing names manually
setnames(pheno.original,c("Movement-WiggleIdex-ExpandedIn'Units'Note-_-100ppm-WI_is_a_ratio:_first_environment_standard_ove_exposed_to_50:50_mix_12C6:13C6_imidacloprid_-99%_analytical_reagent-_concentration_100_ppm_one_hour_Larvae-SexMixed" ,
                          "Movement-WiggleIdex-ExpandedIn'Units'Note-_-25ppm-WI_is_a_ratio:_first_environment_standard_ove_exposed_to_50:50_mix_12C6:13C6_imidacloprid_-99%_analytical_reagent-_concentration_25_ppm_one_hour_Larvae-SexMixed" ) ,
         
         c("Movement-WiggleIdex-ExpandedIn'Units'Note -100ppm-",
           "Movement-WiggleIdex-ExpandedIn'Units'Note-25ppm-")
) 

#make sure to save with these new names
saveRDS(pheno.original, "updatedphenos")
#create vector of phenotype names, exluding first row of lineids
phenos = colnames(pheno.original)[2:length(colnames(pheno.original))]

#save phenotype names
saveRDS(phenos, "updatedphenonames")
phenos = readRDS("updatedphenonames")
#this loop creates directories named after each of the phenotypes, within each folder
foreach(p = unique(phenos))%do%{
  #originalgrms.dir = paste0(ORIGINALOUTPUTGRMS,p,".original.grms")
  permsgrms.dir = paste0(PERMSOUTPUTGRMS,p,".perms.grms")
  
  #originalnogrms.dir = paste0(ORIGINALOUTPUTNOGRMS,p,".original.nogrms")
  permsnogrms.dir = paste0(PERMSOUTPUTNOGRMS,p,".perms.nogrms")
  
  #dir.create(originalgrms.dir)
  dir.create(permsgrms.dir)
  #dir.create(originalnogrms.dir)
  dir.create(permsnogrms.dir)
}



### create gwas id table
#"/project/berglandlab/yangyu/Yang_Adam/GWAS_Input"
#create table that has as colums: phenotype names, permutation #, and task id number
#this will be used during array jobs to decide permutation and phenotype
gwas.table = data.table(expand.grid(pheno = phenos, perms = c(1:10)))
gwas.table[,array.id := c(1:dim(gwas.table)[1])]



write.table(gwas.table,file="/scratch/bal7cg/Yang_Adam/perm.gwas.ref.table",quote=F,row.names=F,col.names=T,sep="\t")
###create wideform tables with permutated line id, and true combination of 
# phenotype line averages- filtered down to our chosen phenotype list
#loading in our two data tables of phenotype data
oldphenos = readRDS("realphenotable")
colnames(oldphenos)[1]= "ral_id"
newphenos = readRDS("./updatesphenos")
#need to add "RAL_ID_" to line names for new phenos
newphenos$ral_id = paste0("RAL_ID_", newphenos$ral_id)
#bring together all phenotypes
totalphenos = merge(oldphenos, newphenos, by = "ral_id")
#awesome, now filter out only the phenotypes we want
phenonames = readRDS("./updatedphenonames")
cols = c("ral_id", phenonames)
selectphenos = totalphenos[,cols,with=F]
#appending in wolbachia status
wol = fread("/project/berglandlab/Yang_Adam/reference_files/female.nStrain20plus.10pca.wolbachia.txt")
wol = wol[,c("FID","wolbachia"),with=F]

setnames(wol,"FID","ral_id")
wol$ral_id = gsub("line", "RAL_ID", wol$ral_id)
selectphenos = merge(selectphenos, wol, by = "ral_id")
#beautiful! save the "real" full phenotable
saveRDS(selectphenos, "combinedphenotable")
#next step- create a set of 10 permutated tables. can create more if more permutations are desired
out = foreach(i = c(1:10)) %do% {
  
  set.seed(i)
  shufflephenos = selectphenos
  shufflephenos$ral_id <- shufflephenos$ral_id[sample(nrow(shufflephenos))]
  write.table(shufflephenos,paste0("./permutated_pheno_tables/permutated_table_", i, ".txt"),quote=F,row.names=F,col.names=T,sep="\t" )
  
}

