### make sliding window figure for joaquin

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
#library(coloc)
#install.packages('doParallel', repos='http://cran.us.r-project.org')
#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(10)
###################################
#gather in the sliding window files
#######################################
### import
fl <- list.files("/scratch/bal7cg/June_2022/objects/sliding.window.outputs/Maxtemp.2.5/", pattern = "out", full.names=T)

gwas.win.o <- foreach(fl.i=fl)%dopar%{
   # fl.i <- fl[1]
   message(fl.i)
   load(fl.i)
   gwas.win.o[,glm.perm.n:=tstrsplit(glm.perm, "/")%>%last]
   gwas.win.o[,glm.perm.n:=tstrsplit(glm.perm.n, "\\.")[[1]]]
   gwas.win.o$glm.perm.n =  sub('.', '', gwas.win.o$glm.perm.n)
   gwas.win.o[,pa:=p.adjust(gwas.win.o$fet.p)]
   
   return(gwas.win.o)
}
gwas.win.o <- rbindlist(gwas.win.o)
saveRDS(gwas.win.o, "/project/berglandlab/Adam/Junemodels.slidingwindow.RDS")



#load data 
dt = readRDS("/project/berglandlab/Adam/Junemodels.slidingwindow.RDS")
#

 dt.grouped = dt %>% 
   group_by(realglm) %>% 
   mutate(p.adjusted = p.adjust(fet.p)) %>% 
   filter(p.adjusted < 0.05) %>%
   as.data.table(.)

 #calculate mean and bounds for permutation data

 dt.window= dt.grouped %>% 
   mutate(perm.st = case_when(realglm == 0 ~ "observed" ,
                              TRUE ~ "permutation")) %>% 
   group_by(start,end,win.i,perm.st, realglm) %>% 
   summarise(N = n()) %>%
   summarise(avg.N = mean(N),
             lowerbound = quantile(N, 0.025),
             upperbound = quantile(N, 0.975)) %>%

   as.data.table(.)
 

 dt.window %>%
   dcast(start + end + win.i ~ perm.st, value.var = "upperbound") %>%
   filter(observed > permutation) %>%
   .$win.i -> sig.wins

 
 final.windows.pos = 
   data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
              mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
   ) %>%
   mutate(start = (mid-0.2)*1e6 ,
          end  = (mid+0.2)*1e6  )
 
#outlier_haplowins = 
#  data.table(win.name = c("inv.start", "win_3.1" ,"win_5.1", "win_9.5" , "inv.stop"),
#             start = c(2125744, 3000026, 5050026, 9500026, 13054180 ),
#             end = c( 2325744, 3150026, 5250026, 9650026, 13254180))

 
 library(magrittr)
 dt.window %<>%
   mutate(sig.v.per = case_when(
     win.i %in% sig.wins ~ "sig",
     TRUE ~ "not.sig"
   ) )
 
 saveRDS(dt.window, "sliding.windowholm.RDS")
 
 
dt.window =  readRDS("sliding.windowholm.RDS")
 ggplot() +
#geom_rect(data=outlier_haplowins, aes(xmin=start/1e6, xmax=end/1e6, ymin=-1, ymax=100),fill ="lightgoldenrod1", alpha=.5) +
   geom_rect(data = final.windows.pos,
             aes(xmin=start/1e6, xmax = end/1e6,
                 ymin = 0, ymax = 40), 
             alpha = 0.7, fill = "gold") +
   #geom_line(data = dt.window[perm.st == "permutation"],
   #           aes(x=I(start/2+end/2)/1e6, y=avg.N), color="black") +
   #geom_line(data = dt.window[perm.st == "permutation"],
   #           aes(x=I(start/2+end/2)/1e6, y=lowerbound), color="black", linetype = "dashed") +
   geom_vline(xintercept = 2225744/1e6) +
   geom_vline(xintercept = 13154180/1e6) +
   geom_ribbon(data = dt.window[perm.st == "permutation"],
             aes(x=I(start/2+end/2)/1e6, ymax= upperbound, ymin = 0), fill="grey", linetype = "solid", alpha = 0.5) +
   geom_point(data=dt.window[perm.st == "observed"],
              aes(x=I(start/2+end/2)/1e6, y=avg.N, color = sig.v.per), #color="firebrick4",
              size = 2) +
   
   theme_bw()+
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos")

 ggsave(g1, file =  "/scratch/bal7cg/new. model2.5.fdr.pdf", w = 7, h = 3)
 
 #we're interested in which phenotypes pass fdr and permutations for each peak - windows 
 dt.window[sig.v.per == "sig"]
 wins.interest = data.table(
   win.i = c(38,39, 88, 97:99)
   #peaks = c("inv.start", "win_3", "win_4.7","win_4.7", "win_5.1", "win_6.2", "win_6.8")
 )
 #find phenotypes that pass fdr correction and are within these windows
 
 dt.grouped = dt %>% 
   group_by(realglm) %>% 
   mutate(p.adjusted = p.adjust(fet.p)) %>% 
   filter(p.adjusted < 0.05) %>%
   as.data.table(.)
 pass.dt = merge(dt.grouped, wins.interest, by = "win.i")
 pass.dt = pass.dt[realglm == 0]
 #we only need window, peak, and phenotype
 pass.dt = pass.dt %>% 
   select( gwas.pheno) %>% 
    distinct(.)
 #save data to use in pca plot