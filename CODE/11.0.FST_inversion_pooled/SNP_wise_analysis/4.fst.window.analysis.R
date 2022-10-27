#### analyse windows of FST
#### 
#### 

### load modules
rm(list = ls())

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)


### data input
args = commandArgs(trailingOnly=TRUE)
analyses.set=as.character(args[1])


#sets=
#"time.eu.w"
#"time.eu.e"
#"time.cville"
#"space.NoA.E"
#"space.eu.w"
#"space.eu.e"
#"core20.seasonality"

######
######
######
######

ag.folder <- "/project/berglandlab/DEST_fst_out/fst.ag.folder"

files.vec = system(paste("ls ",  ag.folder),
                   intern= T)


######## ------> Build function 


#analyses.set="time.eu.e"
#win.bp=1e+05
#step.bp=50000

calculate.win.fst = function(analyses.set, win.bp=10000, step.bp=5000){
  
  
  message(analyses.set)
  
  message("extract types of analyses")
  files.tmp = files.vec[grep(analyses.set , files.vec)]
  
  message(paste(files.tmp, sep = " | "))
  
  tmp.ag = foreach(i = 1:length(files.tmp),
                          .combine = "rbind")%do%{
                            file.tmp = files.tmp[i]
                            message(paste(i, file.tmp, sep = " | "))
                            
                            tmp.load = get(load(paste(ag.folder, file.tmp, sep = "/")))
                            
                            return(tmp.load)
                            
                          }
  
  setDT(tmp.ag)
  
  ## prepare windows
  setkey(tmp.ag, "chr")
  
  message("create window vector")
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L"
                          ,"2R", "3L", "3R"
                          ),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- tmp.ag %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  
  wins[,i:=1:dim(wins)[1]]
  
  dim(wins)
  
  #### --> aggregate
  setkey(tmp.ag, chr, pos)
  head(tmp.ag)
  
  message("aggregate by window -- analysis")
  
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind")%do%{
                       
                       message(paste(win.i, dim(wins)[1], sep=" / "))
                       
                       
                       win.tmp <- tmp.ag[J(data.table(chr=wins[win.i]$chr, 
                                                             pos=wins[win.i]$start:wins[win.i]$end, 
                                                             key="chr,pos")), nomatch=0]
                       
                       win.tmp %>%
                         group_by(comp.set, year_diff, analysis.type) %>%
                         summarise(mean.means.fst.w = mean(mean.SNPwise.FST, na.rm = T),
                                   mean.medians.fst.w = mean(median.SNPwise.FST, na.rm = T)) %>%
                         mutate(chr=wins[win.i]$chr,
                                start=wins[win.i]$start,
                                end=wins[win.i]$end) %>%
                         mutate(mid.point = (start+end)/2 ) ->
                         win.ag.f
                       
                       return(win.ag.f)
                       
                     }
  
  return(win.out)
  message("done")
  
}


######
### analyses.set="core20.seasonality"
tmp.fst.ag = calculate.win.fst(analyses.set, win.bp=2e+05, step.bp=1e+05)

#### plot
#tmp.fst.ag %>%
#  filter(year_diff %in% 0:1) %>%
#  ggplot(aes(x=mid.point,
#             y=mean.means.fst.w,
#             color = as.factor(year_diff)
#             )) +
#  geom_line() +
#  ylim(-0.005,0.005) +
#  geom_hline(yintercept = 0) +
#  facet_grid(as.factor(year_diff)~chr, scales = "free_x") ->
#  fst.ag.plot
#
#ggsave(fst.ag.plot,
#       file =  paste(analyses.set, "fst.plot", "pdf", sep = "."),
#       h= 3 , w = 9)

save(tmp.fst.ag,
     file = paste(analyses.set, "window.fst.ag", "Rdata", sep = "."))

