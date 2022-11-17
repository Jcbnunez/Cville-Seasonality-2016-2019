###aggregate data into 1 minute sumes
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
#library(ggrepel)
registerDoParallel(5)
#rivanna wd
#wd = "/home/bal7cg/R/data.objects"
wd = "/scratch/bal7cg/Inverson_project/"
wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
setwd( wd )

week1 = readRDS("week1.def.data")
week2 = readRDS("week2.def.data")
week3 = readRDS("week3.def.data")
#add week to fly.id
week1$fly.id = paste0(week1$fly.id,"_week1" )
week2$fly.id = paste0(week2$fly.id,"_week2" )
week3$fly.id = paste0(week3$fly.id,"_week3" )
total.def.data = rbind(week1, week2, week3)
saveRDS(total.def.data, "total.def.data")

total.def.data = readRDS("total.def.data")
head(total.def.data)



#ditch unneeded columns for ease
total.def.data = total.def.data %>% 
  dplyr::select( !c(fly.geno, wells, monitor, well.number, day))

# #how many flies have atleast 15 minutes post stim data
# topdata = total.def.data[cmin >= 20]
# length(unique(topdata$fly.id))
#create a vector of flyids
fly.ids = total.def.data$fly.id #yeah i don't know, this wasn't working all on the same line
fly.ids = unique(fly.ids)
#create a table that defines a set of windows from -62 to 30, window size 90 seconds, 30 second step. 
startpoints = c(seq(-62, -3.5, 0.5), seq(0, 28.5, 0.5))
endpoints =  c( seq(-60.5,-2, 0.5), seq(1.5, 30, 0.5))

flyid.out = foreach(f = fly.ids, .packages = c("data.table", "tidyverse","foreach")) %dopar% {
  f = fly.ids[1406]
  #filter out data
  fly.data = total.def.data %>% 
    filter(fly.id == f)
  b.a = unique(fly.data$basal_act)
  #basal activity here is basal_act / 60 , for activity level over 1 minute
  
  sliding.window.o = foreach(g = 1:length(startpoints), .combine = "rbind") %do% {
    #g = 76
    start.pt = startpoints[g]
    end.pt = endpoints[g]
    #filter to only include time points within these time points
    time.data = fly.data[cmin >= start.pt & cmin <= end.pt]
    window.act = sum(time.data$activity)
    unique(time.data$basal_act)
    #output window #, stpoitn, endpoint, and sum activity
    #check to see if this time frame has or is missing data- for missing data report na
    if (dim(time.data)[1] == 0) {
      o = data.frame(
        minute = (start.pt + end.pt) / 2,
        start.point = start.pt,
        end.point = end.pt,
        window.activity = NA
      )
    } else {
      o = data.frame(
        minute = (start.pt + end.pt) / 2,
        start.point = start.pt,
        end.point = end.pt,
        window.activity = window.act
      )
    }
    
    #find first time point below b.a
    o
  } 
  # ggplot(sliding.window.o, aes(x = minute,  y = window.activity)) + geom_line() +
  #   geom_hline(yintercept = (b.a / 60)) +
  # theme_bw() +
  #   xlim(-20,20)
  
  sliding.window.o$fly.id = f
  return(sliding.window.o)
  
  
}
fly.long = rbindlist(flyid.out)
saveRDS(fly.long, "data.smoothed90window30step.RDS")

