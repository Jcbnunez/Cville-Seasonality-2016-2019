#### Omnibus enrrichment for windows
#### 

#args = commandArgs(trailingOnly=TRUE)
#job=as.numeric(args[1])-1
#k=as.numeric(args[1])-1
#model=args[2]
#k=0


### libraries
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)

####

####  create window objects
final.windows.pos = 
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

### load thermal GLM object
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

models = c(
		   #"temp.max;4;5.Cville",
		   "humidity.ave;10;1.Europe_W",
		   "humidity.var;3;3.Europe_E",
		   #"temp.min;8;2.North_America_E",
		   "temp.min;8;2.North_America_E."
		    #,
           #"temp.ave;9;3.Europe_E" #,
           #"temp.ave;1;2.North_America_E" #,
           #"humidity.ave;8;1.Europe_W"
           )

base <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper"


######
######
### ----> loac Cville object ... this is a constant througt
#file.cvile <- paste(base, "temp.max;2;5.Cville", "temp.max;2;5.Cville.glmRNP.Rdata", sep = "/" )
#print(file.cvile)
file.cvile <- paste(base, "Revision_Best_Models", paste( "temp.max;2;5.Cville","v2.glmRNP.Rdata", sep = ".") , sep = "/" )
  print(file.cvile)

out.glm.cvile <- get(load(file.cvile))

enrichment.sets =
foreach(z = 1:length(models)
        ,.combine = "rbind"
        )%do%{
 
  model = models[z]
  message(model)
  
  ### ----> load relative object
base <- "/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper"

  file <- paste(base, "Revision_Best_Models", paste(model,"v2.glmRNP.Rdata", sep = ".") , sep = "/" )
  
  print(file)
  out.glm.matched <- get(load(file)) 
  
  k = 0
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
  m1[!is.na(p_lrt.x.x) & !is.na(p_lrt.x.y) ,cville.rank := rank(p_lrt.x.x)/length(p_lrt.x.x)]
  m1[!is.na(p_lrt.x.x) & !is.na(p_lrt.x.y) ,anchor.rank := rank(p_lrt.x.y)/length(p_lrt.x.y)]
  
  message("Estimating the betas")
  
  m1.beta <- m1[,list(cville.beta=b_temp.x,
                      anchor.beta=b_temp.y), 
                list(variant.id)]
  
  m1 <- merge(m1, m1.beta, by="variant.id")
  m1[,inv:=invName.y!="none"]
  
  thrs <- c(
    0.05
  )
  
  wins=final.windows.pos
  
  o.win <- foreach(#chr.i=unique(m1$chr.x), 
                    chr.i="2L",
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
                                         
                                         ###  enrichment
                                         tab <- table(m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$cville.rank  < thr.i,
                                                      m1[chr.x==chr.i][ pos.x>=wins$start[pos.ith] & pos.x<=wins$end[pos.ith] ]$anchor.rank <  thr.i)
                                         
                                         fet <- fisher.test(tab)
                                         
                                         ###  sign test
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
  
  #o.win %>%
  #  filter(p < 0.05)


} ## close model loop

message("Save File")

##best_model =
##c("temp.ave;9;3.Europe_E",
##"humidity.ave;8;1.Europe_W",
##"temp.ave;1;2.North_America_I95"#,
###"humidity.ave;11;2.North_America_Mid"
##)

##temp_model =
##c("temp.max;2;3.Europe_E",
##"temp.max;2;1.Europe_W",
##"temp.max;2;2.North_America_I95" #,
###"temp.max;2;2.North_America_Mid"
##)

#enrichment.sets %<>%
#  mutate(analysis_type = case_when(anchor.model %in% best_model ~ "best_model",
#                                   anchor.model %in% temp_model ~ "temp_model"
#  ))

save(enrichment.sets, 
     file = "./window.enrich.set.Rdata"
      # paste("out.enr/",
      #       paste(model, k , "Omnubis.enrich.winlevel", "Rdata",
      #             sep = "."),
      #       sep = "")
)    



    