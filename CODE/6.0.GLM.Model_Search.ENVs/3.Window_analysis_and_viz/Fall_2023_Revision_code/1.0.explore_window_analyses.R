### Analysis of New models

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

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
sets

#####
### open GDS
genofile <- seqOpen("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

###
### samps
samps <- fread("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

###
#mods_fin$mod_var -> models

models = c("temp.max;2;5.Cville",
		   "temp.max;4;5.Cville",
		   "humidity.ave;10;1.Europe_W",
		   "humidity.var;3;3.Europe_E",
		   "temp.min;8;2.North_America_E"
		    #,
           #"temp.ave;9;3.Europe_E" #,
           #"temp.ave;1;2.North_America_E" #,
           #"humidity.ave;8;1.Europe_W"
           )

#### base files
base <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/"

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
#k=1

#all.wins.out = 
#foreach(k=1:length(models), .combine = "rbind")%do%{
#

file <- paste(base, "Revision_Best_Models", paste( models[k],"v2.glmRNP.Rdata", sep = ".") , sep = "/" )
  print(file)
  
  message(models[k])
  out.glm <- get(load(file))
  
  ##############
  ##############
  
  ###############
  ### windows ###
  ###############
  # generate a master index for window analysis
  ### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  
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
  
  ### run windows
  
  setkey(out.glm, chr, pos)
  head(out.glm)
  
  #wins %<>% filter(chr == "2L")
  
  ### start the summarization process
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind"
  )%do%{
    
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.tmp <- glm.out[J(data.table(chr=wins[win.i]$chr, 
                                    pos=wins[win.i]$start:wins[win.i]$end, 
                                    key="chr,pos")), nomatch=0]
    
    
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt.x, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]
    
    
    seqSetFilter(genofile, 
                 variant.id=unique(win.tmp$variant.id),
                 sample.id=samps[locality=="VA_ch"][year>= 2016 ]$sampleId)
    
    #obtain AFs 
    af <- seqGetData(genofile, "annotation/format/FREQ")
    f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), 
                        variant.id=seqGetData(genofile, "variant.id"))
    
    #merge AFs with object
    win.tmp <- merge(win.tmp, f.hat, by="variant.id")
    win.tmp[,het:=2*fhat*(1-fhat)]
    
    #thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))[,1]
    
    pr.i <- c(
              0.05
              )
    
    #tmpo <- foreach(
    # pr.i=thrs, 
    # .errorhandling="remove", 
    # .combine="rbind")%do%{
    win.tmp %>% 
      filter(!is.na(rnp), pr.i == pr.i ) %>%
      group_by(perm, chr , variable, mod) %>%
      summarise(pos_mean = mean(pos),
                pos_mean = mean(pos),
                pos_min = min(pos),
                pos_max = max(pos),
                win=win.i,
                pr=pr.i,
                rnp.pr=c(mean(rnp<=pr.i)),
                rnp.binom.p=c(binom.test(sum(rnp<=pr.i), 
                                         length(rnp), pr.i)$p.value),
                wZa=sum(het*Z)/(sqrt(sum(het^2))),
                wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                min.p.lrt=min(p_lrt.x),
                min.rnp=min(rnp),
                nSNPs = n(),
                sum.rnp=sum(rnp<=pr.i),
      ) %>%
      mutate(
        model.pop = models[k],
        perm_type=ifelse(perm==0, "real","permuted"),
        invName=case_when(
          chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
          chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
          chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
          chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
          TRUE ~ "noInv"
        )) -> win.out
        
    #}
    #tmpo
  }

##}

### save
out_folder <- "/gpfs2/scratch/jcnunez/genetics_resub/7.newmodels"

message(paste(out_folder, "/Window_analysis_", models[k] , ".Rdata", sep=""))
save(win.out, 
     file=paste(out_folder, "/Window_analysis_", models[k] , ".Rdata", sep="")
     )
