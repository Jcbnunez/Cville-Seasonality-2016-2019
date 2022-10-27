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
library(SeqArray)
library(doMC)
library(car)
library(DescTools)
library(ape)
library(ggtree)
library(aplot)
library(forcats)
registerDoMC(2)

vcf_file <- "./Spn27A.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )
x <- vcfR2genlight(vcf)

tab(x) %>%
  as.data.frame() %>%
  filter(`2L_6671911_SNP` %in% c(0,1,2) ) %>%
  mutate(sampleId = rownames(.)) %>%
  .[grep("OW", .$sampleId ),] %>%
  separate(sampleId, into = c("Exp", "Well", "PopSeason"), sep = "_") ->
  marker.dat

marker.dat %>%
  ggplot(aes(x= PopSeason, 
             fill= as.factor(`2L_6671911_SNP`))) +
  geom_bar(position = "dodge") ->
  bar_slc
ggsave(bar_slc, file ="bar_slc.pdf")

marker.dat %>%
  group_by(PopSeason, `2L_6671911_SNP`) %>%
  summarise(N = n())

Input =("
 Allele  Spring  Fall      
 0         26      31
 1         8        1
")

Mat = as.matrix(read.table(textConnection(Input),
                              header=TRUE,
                              row.names=1))

Mat

library(DescTools)
GTest(Mat)

