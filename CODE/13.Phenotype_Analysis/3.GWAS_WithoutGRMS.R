### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(ggplot2)














args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = args[1]

id = 1 #use this when testing to see if script functions
### read in data files

INPUT = "/scratch/bal7cg/Yang_Adam/"
ingds = paste0(INPUT,"dgrp.gds")
#new.gwas.ref.table.txt is the task id file for new phenotypes
#gwas.ref.table.txt is the task id file for old phenotypes
#updatedphenos is the phenotype data for new phenotypes
#realphenotable is the phenotype data for old phenotypes


task.id = fread(paste0(INPUT,"new.gwas.ref.table.txt"))#load in task id for observed data
perms.tables = readRDS(paste0(INPUT,"updatedphenos"))

#need to ensure the lineid column matches with gds file
perms.tables$ral_id <- paste("line_", perms.tables$ral_id, sep = "" )

phenotypes = task.id[array.id == id]$pheno #define phenotype based on task id


cols = c("ral_id",phenotypes)

dt = perms.tables[,cols,with=F]
#load in wolbachia status
wol = fread("/project/berglandlab/Yang_Adam/reference_files/female.nStrain20plus.10pca.wolbachia.txt")
wol = wol[,c("FID","wolbachia"),with=F]

setnames(wol,"FID","ral_id")

setkey(dt,ral_id)
setkey(wol,ral_id)
#bind phenotype data with wolbachia status
dt = merge(dt,wol)




#re organize
colnames(dt)[2] <- "avg"
#create an "identity matrix" to represent the GRM, has dimensions equal to number of lines, and row and column names matching the line ids
GRM = diag(dim(dt[!is.na(avg)])[1])
colnames(GRM) <- dt[!is.na(avg)]$ral_id
rownames(GRM) <- dt[!is.na(avg)]$ral_id



#fit a GLMM to the data
#GMMAT creates a null model that uses wolbachia status as fixed effect on phenotype, and the identity matrix as a (nonexistent) random effect. family refers to the method used to calculate the model
modelqtl <- glmmkin(fixed =avg ~ wolbachia, data = dt[!is.na(avg)], kins = GRM, id = "ral_id",family = gaussian(link = "identity"))



#run GMMAT
#create  path for output gmmat file to be saved to- using the phenotype variable allows it to be saved to the corresponding directory
outputs = paste("/scratch/bal7cg/Yang_Adam/GWAS_withoutGRMs/NewPhenoOutputsWithoutGRMS/" , 
                phenotypes, ".original.nogrms/", 
                phenotypes,".nogrms.original.txt",
                sep = "")


#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.5), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 4)






print ("done")


