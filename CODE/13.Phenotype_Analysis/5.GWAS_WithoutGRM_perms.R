### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(ggplot2)




args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = args[1]

#id = 1

INPUT = "/scratch/bal7cg/Yang_Adam/"
ingds = paste0(INPUT,"dgrp.gds")
#load in the task.id table, use it to define permutation and phenotype
task.id = fread(paste0(INPUT,"perm.gwas.ref.table.txt"))
phenotypes = task.id[array.id == id]$pheno
permutations = task.id[array.id == id]$perms


#load in the permutation table

perms.tables = fread(paste0(INPUT,"permutated_pheno_tables/permutated_table_",permutations,".txt"))
#define columns we want to keep
cols = c("ral_id",phenotypes,"wolbachia")


dt = perms.tables[,cols,with=F]



#re organize
colnames(dt)[2] <- "avg"
#have to change the id names to match grm
dt$ral_id = gsub("RAL_ID", "line", dt$ral_id)


#create grm as an identity matrix
grm = diag(dim(dt[!is.na(avg)])[1])
colnames(grm) <- dt[!is.na(avg)]$ral_id
rownames(grm) <- dt[!is.na(avg)]$ral_id

#fit a GLMM to the data
modelqtl <- glmmkin(fixed =avg ~ wolbachia, data = dt[!is.na(avg)], kins = grm, id = "ral_id",family = gaussian(link = "identity"))



#run GMMAT
outputs = paste("/scratch/bal7cg/Yang_Adam/GWAS_withoutGRMs/PermsOutputsWithoutGRMs/" , 
                phenotypes, ".perms.nogrms/", 
                phenotypes,".nogrms.perms." , permutations, ".txt",
                sep = "")



glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, 
           nperbatch = 400, ncores = 4)





print ("done")

