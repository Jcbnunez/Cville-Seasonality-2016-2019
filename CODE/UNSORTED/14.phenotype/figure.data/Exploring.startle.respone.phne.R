### 
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(patchwork)

GWAS.sr <- fread("/project/berglandlab/Yang_Adam/YangsGWAS/GWAS_withoutGRMs/Cynthiaphenotypes/StartleResponse_standard_Female.original.nogrms/StartleResponse_standard_Female.nogrms.original.txt")



GWAS.sr %>%
  filter(CHR == "2L") ->
  GWAS.sr.2l

GWAS.sr.2l %>%
  filter(PVAL < 0.01) %>%
  arrange(-abs(SCORE)) %>%
  ggplot(
    aes(
      x=POS,
      y=-log10(PVAL)
    )
  ) + geom_point() +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  signif.beta

ggsave(signif.beta, file = "signif.beta.pdf", w = 6, h = 3)
  

###  ### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(GWAS.sr.2l, "CHR")

## prepare windows
wins <- foreach(chr.i=c("2L"
                        #,"2R", "3L", "3R"
                        ),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- GWAS.sr.2l %>%
                    filter(CHR == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$POS), to=max(tmp$POS)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
#####
#####
#####

win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- GWAS.sr.2l %>%
    filter(CHR == wins[win.i]$chr) %>%
    filter(POS >= wins[win.i]$start & POS <= wins[win.i]$end)
  
  data.frame(
    CHR =wins[win.i]$chr,
    START =wins[win.i]$start, 
    END = wins[win.i]$end,
    MEAN_VAR = mean(win.tmp$VAR),
    MEAN_ABS_SCORE = mean(abs(win.tmp$SCORE)),
    MEAN_RAW_SCORE = mean((win.tmp$SCORE)),
    MEAN_PVAL = mean(-log10(win.tmp$PVAL)),
    NUM_Pval_0.05 = dim(filter(win.tmp,  PVAL < 0.05 ))[1],
    NUM_Pval_0.01 = dim(filter(win.tmp,  PVAL < 0.01 ))[1]
    #NUM_Pval_0.0001 = dim(filter(win.tmp,  PVAL < 0.0001 ))[1]
  )
}

win.out %>%
  reshape2::melt(id = c("CHR",  "START",    "END")) %>% 
  ggplot(
    aes(
      x=c(START+END)/2,
      y=value
    )
  ) + geom_line() +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) +
  facet_wrap(~variable, scales = "free") ->
  startle.plot.mean

ggsave(startle.plot.mean, file ="startle.plot.mean.pdf", w = 9, h = 6)



