### Draw haplonets
### 

library(pegas)
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)
library(foreach)
library(doMC)
registerDoMC(2)

setwd("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/")

# load outlier GLMs
outliers <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

#####
load("./AF_dat_summ_id.Rdata")
#---> AF_dat_summ_id
win5 = filter(AF_dat_summ_id, win == "win_5", pos %in% outliers$pos )
win6 = filter(AF_dat_summ_id, win == "win_6", pos %in% outliers$pos)
win10 = filter(AF_dat_summ_id, win == "win_10", pos %in% outliers$pos)

####

STD_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Std_samps_OnlyNames.txt", header = F)
STD_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/STD_DGRP_OnlyNames.txt", header = F)

INV_CM = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/Inv_samps_OnlyNames.txt", header = F)
INV_DGRP = fread("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/INV_DGRP_OnlyNames.txt", header = F)

inv_samps = c(
paste(c(INV_CM$V1, INV_DGRP$V1), "_0", sep = ""),
paste(c(INV_CM$V1, INV_DGRP$V1), "_1", sep = ""))

inv_samps %<>%
  data.frame(hap = .) %>%
  mutate(pop=   case_when(hap %in%  paste(c(INV_CM$V1), "_0", sep = "") ~ "CM",
                     hap %in%  paste(c(INV_CM$V1), "_1", sep = "") ~ "CM",
                     hap %in%  paste(c(INV_DGRP$V1), "_0", sep = "") ~ "DGRP",
                     hap %in%  paste(c(INV_DGRP$V1), "_1", sep = "") ~ "DGRP" ))

vcf_file <- "./MC.DGRP.merged.forHap.maf.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE)

inv_ids = which(rownames(my_dnabin1) %in% inv_samps$hap)

bin_win5 = my_dnabin1[inv_ids,which(colnames(my_dnabin1) %in% win5$SNP_id)]
bin_win6 = my_dnabin1[inv_ids,which(colnames(my_dnabin1) %in% win6$SNP_id)]
bin_win10 = my_dnabin1[inv_ids,which(colnames(my_dnabin1) %in% win10$SNP_id)]

#####

h <- haplotype(bin_win5, strict = TRUE)
sz <- summary(h)
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = TRUE)
nt.labs <- attr(nt, "labels")
sz <- sz[nt.labs]

R <- haploFreq(bin_win5, fac = inv_samps$pop, haplo = h)
R <- R[nt.labs, ]

pdf("hap.win5.pdf")
plot(nt, size = sz, pie = R,  labels = TRUE, fast = FALSE, threshold = 0, scale.ratio = 0.9, cex = 0.6, show.mutation = 2 )
dev.off()

#pdf("hap.win5.muts.pdf")
#plot(nt, size = sz, pie = R, show.mutation = 3, labels = FALSE, legend = c(-100, 100))
#dev.off()

######


h <- haplotype(bin_win6)
sz <- summary(h)
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = TRUE)
nt.labs <- attr(nt, "labels")
sz <- sz[nt.labs]

R <- haploFreq(bin_win5, fac = inv_samps$pop, haplo = h)
R <- R[nt.labs, ]

pdf("hap.win6.pdf")
plot(nt, size = sz, pie = R,  labels = TRUE, fast = FALSE, threshold = 0, scale.ratio = 0.8, cex = 0.6, show.mutation = 2 )
dev.off()

####
####


h <- haplotype(bin_win10)
sz <- summary(h)
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = TRUE)
nt.labs <- attr(nt, "labels")
sz <- sz[nt.labs]

R <- haploFreq(bin_win5, fac = inv_samps$pop, haplo = h)
R <- R[nt.labs, ]

pdf("hap.win10.pdf")
plot(nt, size = sz, pie = R,  labels = TRUE, fast = FALSE, threshold = 0, scale.ratio = 0.7, cex = 0.6, show.mutation = 2 )
dev.off()


#### draw joint 
#### 
#### 
#### 

vcf_file <- "./MC.DGRP.merged.forHap.maf.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE)

inv_ids = which(rownames(my_dnabin1) %in% inv_samps$hap)


bin_joint =  my_dnabin1[,which(colnames(my_dnabin1) %in% c(win6$SNP_id,win10$SNP_id) )]


h <- haplotype(bin_joint)
sz <- summary(h)
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = F)
nt.labs <- attr(nt, "labels")
sz <- sz[nt.labs]

R <- haploFreq(bin_win5, fac = inv_samps$pop, haplo = h)
R <- R[nt.labs, ]

pdf("hap.join.6n10.pdf")
plot(nt, size = sz, pie = R,  labels = TRUE, fast = FALSE, threshold = 0, scale.ratio = 0.9, cex = 0.6, show.mutation = 2 )
dev.off()


