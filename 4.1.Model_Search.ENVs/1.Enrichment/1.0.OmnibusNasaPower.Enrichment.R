# ijob -A berglandlab_standard -c1 -p standard --mem=8G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

args = commandArgs(trailingOnly=TRUE)
#job=as.numeric(args[1])-1
k=as.numeric(args[1])-1
model=args[2]

k=0
model="temp.ave;9;3.Europe_E"

### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)

### load thermal GLM object
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

### load core20 GLM object
message("load other pops  GLM")
#load(file="/project/berglandlab/alan/core20glm.Rdata")


base <- "/project/berglandlab/alan/environmental_ombibus_global"
file <- paste(base, model, paste(model,"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)
out.glm.matched <- get(load(file))

###
file.cvile <- paste(base, "temp.max;2;5.Cville", "temp.max;2;5.Cville.glmRNP.Rdata", sep = "/" )
print(file.cvile)
out.glm.cvile <- get(load(file.cvile))

####
out.glm.matched %>%
  filter(perm == k) -> out.glm.matched.mod

out.glm.cvile %>%
  filter(perm == k) -> out.glm.cvile.mod


### merge
message("Merge datasets")

m1 <- merge(out.glm.cvile.mod, 
            out.glm.matched.mod, 
            by.x="variant.id", by.y="variant.id")

m1 <- m1[chr.x!="X"]
m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,cville.rank := rank(p_lrt.x)/length(p_lrt.x)]
m1[!is.na(p_lrt.x) & !is.na(p_lrt.y) ,anchor.rank := rank(p_lrt.y)/length(p_lrt.y)]

message("Estimating the betas")

m1.beta <- m1[,list(cville.beta=b_temp.x,
                    anchor.beta=b_temp.y), 
              list(variant.id)]

m1 <- merge(m1, m1.beta, by="variant.id")
m1[,inv:=invName.y!="none"]

thrs <- c(
  #0.01 #, 
  #0.002
  0.05
  )


message("Running Enrichment")

o <- foreach(chr.i=unique(m1$chr.x), 
             .combine="rbind", 
             .errorhandling="remove")%do%{
               
               foreach(inv.i=c(T,F), 
                       .combine="rbind", 
                       .errorhandling="remove")%do%{
                         
                         foreach(thr.i=thrs, 
                                 .combine="rbind", 
                                 .errorhandling="remove")%do%{
                                   
                                   # chr.i <- "3R"; inv.i <- T; thr.i<-0.05
                                   message(paste(chr.i, inv.i, thr.i, sep=" / "))
                                   
                                   ### Machado set, enrichment
                                   tab <- table(m1[chr.x==chr.i][inv==inv.i]$cville.rank  < thr.i,
                                                m1[chr.x==chr.i][inv==inv.i]$anchor.rank <  thr.i)
                                   
                                   fet <- fisher.test(tab)
                                   
                                   ### Machado set, sign test
                                   st.T <- sum(sign(m1[chr.x==chr.i][inv==inv.i][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) ==
                                                 sign(m1[chr.x==chr.i][inv==inv.i][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                                   st.F <- sum(sign(m1[chr.x==chr.i][inv==inv.i][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) !=
                                                 sign(m1[chr.x==chr.i][inv==inv.i][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                                   
                                   bt <- binom.test(st.T, st.T+st.F, .5)
                                   
                                   
                                   tmp1 <- data.table(chr.x=chr.i, inv=inv.i, thr=thr.i , 
                                                      perm=k,
                                                      anchor.model=model,
                                                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                                                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])
                                   
                                   
                                   
                                   ### return
                                   return(tmp1)
                                   
                                 }
                       }
             }
o

message("Save File")

save(o, 
     file = 
       paste("out.enr/",
       paste(model, k,  "Omnubis.enrich.CHRlevel", "Rdata",
             sep = "."),
       sep = ""
       )
)    

##### Window level analysis
### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(out.glm.cvile.mod, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- out.glm.cvile.mod %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)

#wins %<>% filter(chr == "2L")
o.win <- foreach(chr.i=unique(m1$chr.x), 
              #chr.i="2L",
             .combine="rbind", 
             .errorhandling="remove")%do%{
               
               foreach(pos.ith=1:dim(filter(wins, chr == chr.i))[1], 
                       .combine="rbind", 
                       .errorhandling="remove")%do%{
                         
                         foreach(thr.i=thrs, 
                                 .combine="rbind", 
                                 .errorhandling="remove")%do%{
                                   
                                   # chr.i <- "3R"; inv.i <- T; thr.i<-0.05
                                   message(paste(chr.i, pos.ith, dim(filter(wins, chr == chr.i))[1] , thr.i, sep=" / "))
                                   
                                   ### Machado set, enrichment
                                   tab <- table(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$cville.rank  < thr.i,
                                                m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$anchor.rank <  thr.i)
                                   
                                   fet <- fisher.test(tab)
                                   
                                   ### Machado set, sign test
                                   st.T <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) ==
                                                 sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                                   st.F <- sum(sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$cville.beta) !=
                                                 sign(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ][cville.rank<thr.i & anchor.rank<thr.i]$anchor.beta))
                                   
                                   bt <- binom.test(st.T, st.T+st.F, .5)
                                   
                                   
                                   tmp1 <- data.table(chr.x=chr.i, 
                                                      thr=thr.i , 
                                                      perm=k,
                                                      anchor.model=model,
                                                      win.start=wins$start[pos.ith],
                                                      win.end=wins$end[pos.ith],
                                                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                                                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])
                                   
                                   ### return
                                   return(tmp1)
                                   
                                 }
                       }
             }
o.win %>%
  filter(p < 0.05)


message("Save File")

save(o.win, 
     file = 
       paste("out.enr/",
       paste(model, k , "Omnubis.enrich.winlevel", "Rdata",
             sep = "."),
       sep = "")
)    
