args <- commandArgs(trailingOnly=T)
suffix1 <- as.character(args[2])
prefix <- as.character(args[3])
num <- as.numeric(args[4])
source("/home/nzx3cc/baypass_2.4/utils/baypass_utils.R")
pbc <- read.table(paste0(prefix, "_summary_beta_params.out"),h=T)$Mean
data <- geno2YN(paste0("/scratch/nzx3cc/nzx3cc/rawdata_baypass/", prefix, ".genobaypass"))
if (prefix == "all"){
omega = as.matrix(read.table("/scratch/nzx3cc/nzx3cc/rawdata_baypass/allthinned_mat_omega.out"))
} else{
omega = as.matrix(read.table(paste0("/scratch/nzx3cc/nzx3cc/rawdata_baypass/no", prefix, "thinned_mat_omega.out")))}
simulate.baypass(omega.mat = omega, nsnp = num, sample.size=data$NN, beta.pi=pbc, pi.maf=0, suffix=suffix1)
