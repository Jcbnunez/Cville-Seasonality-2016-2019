library(tidyverse)
library(magrittr)
library(vroom)

library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)
library(patchwork)
library(foreach)
library(doMC)
library(doParallel)
registerDoMC(4)


args <- commandArgs(trailingOnly = TRUE)
file1 <- as.character(args[2])
bayes <- read.table(file1, header=T)
win.bp <- 1e5
step.bp <- 5e4

wins <- 

foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- bayes %>%
                      filter(CHR == chr.i, COVARIABLE==1)
                    
                    data.table(CHR=chr.i,
                               start=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp) + win.bp)
                 }
wins[,i:=1:dim(wins)[1]]

print(dim(wins))
