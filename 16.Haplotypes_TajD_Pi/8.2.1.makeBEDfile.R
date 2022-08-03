### make bed file 
### 
rm(list = ls())

### libraries
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


models = c("temp.max;2;5.Cville" #,
           #"temp.ave;9;3.Europe_E",
           #"humidity.ave;4;2.North_America_W",
           #"humidity.ave;8;1.Europe_W"
)

#### base files
base <- "/project/berglandlab/alan/environmental_ombibus_global"
k=1

#all.wins.out = 
#foreach(k=1:length(models), .combine = "rbind")%do%{
#

file <- paste(base, models[k], paste(models[k],"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)

message(models[k])
out.glm <- get(load(file))


# generate a master index for window analysis
### define windows
win.bp <- 1e5
step.bp <- 1e5+1

setkey(out.glm, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- glm.out %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
#####

## BED file in2l
wins %>% filter(chr == "2L") %>% 
  mutate(start.bed = start-1, end.bed = end-1) %>% 
  dplyr::select(chr, start.bed,  end.bed) ->
  bed.file.2L

write.table(bed.file.2L, 
            file = "bed_guide.file.2L.BED", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")