### libraries
library(data.table)
library(SeqArray)
library(rehh)
library(tidyverse)
library(foreach)
library(car)
library(DescTools)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

setwd("/scratch/yey2sn/Overwintering_ms/13.Selection_REHH")

load("CM.IHS.HH.INV.out.Rdata")

wgscan.ihs$ihs %>%
  as.data.frame() %>%
  #mutate(POSITION_bin = RoundTo(POSITION, 10000, "floor")) %>%
  #group_by(POSITION_bin) %>%
  #summarize(Mean_ihs = mean(abs(IHS))) %>%
  ggplot(aes(
    #x=POSITION_bin,
    #y=Mean_ihs,
    x=POSITION,
    y=(IHS)
    #color = LOGPVALUE
  )) + 
  geom_point() +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  ihs_point

ggsave(ihs_point, file = "ihs_point.png")

  
  
  
