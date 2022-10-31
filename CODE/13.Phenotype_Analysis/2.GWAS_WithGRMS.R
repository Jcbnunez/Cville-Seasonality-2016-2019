### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(ggplot2)



args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = args[1]

#id = 1#use this when testing to see if script functions correctly


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/scratch/bal7cg/Yang_Adam/"#path of working directory
ingds = paste0(INPUT,"dgrp.gds")
#load in the task.id table created for observed (not permutations), use it to define  phenotype
task.id = fread(paste0(INPUT,"new.gwas.ref.table.txt"))#load in task id for observed data
#in this instance, 
phenotypes = task.id[array.id == id]$pheno

#load in the permutation table

perms.tables = readRDS(paste0(INPUT,"updatedphenos"))
#define columns we want to keep
cols = c("ral_id",phenotypes,"wolbachia")
#load in grm (from dgrp website)
grm = readRDS(paste0(INPUT, "GRM.Aug"))

dt = perms.tables[,cols,with=F]

#load in wolbachia status data, merge with phenotype data
wol = fread("/project/berglandlab/Yang_Adam/reference_files/female.nStrain20plus.10pca.wolbachia.txt")
wol = wol[,c("FID","wolbachia"),with=F]

setnames(wol,"FID","ral_id")

setkey(dt,ral_id)
setkey(wol,ral_id)

dt = merge(dt,wol)


#re organize
colnames(dt)[2] <- "avg"
#have to change the id names to match grm
dt$ral_id = gsub("RAL_ID", "line", dt$ral_id)


#GMMAT creates a null model that uses wolbachia status as fixed effect on phenotype, and the identity matrix as a (nonexistent) random effect. family refers to the method used to calculate the model
#fit a GLMM to the data
modelqtl <- glmmkin(fixed =avg ~ wolbachia, data = dt[!is.na(avg)], kins = grm, id = "ral_id",family = gaussian(link = "identity"))



#run GMMAT
#create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
outputs = paste("/scratch/bal7cg/Yang_Adam/GWAS_withGRMs/PermsOutputsWithGRMS/" , 
                phenotypes, ".perms.grms/", 
                phenotypes,".grms.perms." , permutations, ".txt",
                sep = "")


#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 4)


print ("done")


