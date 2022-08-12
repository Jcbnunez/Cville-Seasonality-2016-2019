library(magrittr)
library(tidyverse)
library(vroom)

setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/4.3.Case.Study.Msp300/msp.300.structure/")
### Mps 300 gff
gff.msp.300 = vroom("gff.msp.300.dat.txt", col_names = F, delim = "\t")
for(i in 1:dim(gff.msp.300)[1] ){
  tmp <- gff.msp.300[i,]
  
  tmp.str = str_split(tmp, pattern = "\\|")
  
  isoform = gsub('.+(isoform .).+','\\1',tmp.str[[9]][3])
  
  gff.msp.300$isoform[i] = isoform
  
}
gff.msp.300$isoform %>% unique()

gff.msp.300 %<>% 
  mutate(rank = case_when(isoform == "isoform B" ~ -0.001,
                          isoform == "isoform D" ~ -0.002,
                          isoform == "isoform E" ~ -0.003,
                          isoform == "isoform F" ~ -0.004,
                          isoform == "isoform G" ~ -0.005,
                          isoform == "isoform H" ~ -0.006,
                          isoform == "isoform I" ~ -0.007,
                          isoform == "isoform J" ~ -0.008,
                          isoform == "isoform K" ~ -0.009,
                          isoform == "isoform L" ~ -0.010,
                          isoform == "isoform M" ~ -0.011
  ))


gff.msp.300 %>%
  #filter(isoform == "isoform D") %>%
  mutate(dist= X5-X4) -> exons

exons %<>% 
  mutate(exon.of.int = case_when(X4 < 5192177 & X5 > 5192177 ~ "contains.mut",
                                 X4 > 5204514 & X5 < 5205349 ~ "contains.kash",
                                 TRUE ~ "none"))

exons %>% filter(exon.of.int != "none")

