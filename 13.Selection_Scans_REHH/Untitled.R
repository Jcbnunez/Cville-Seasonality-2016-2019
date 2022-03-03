### libraries
library(data.table)
library(SeqArray)
library(tidyverse)
library(foreach)
library(car)
library(DescTools)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

load("CM.IHS.HH.out.Rdata")


