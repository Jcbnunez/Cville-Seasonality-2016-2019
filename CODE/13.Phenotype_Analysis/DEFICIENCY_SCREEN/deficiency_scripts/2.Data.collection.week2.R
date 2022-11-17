#examine extent of startle response recorded using Trikinetics
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggrepel)

#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(10)
#create a list that contains the unique values for each day- stimulus times, and monitors used, and wd
day1list = list(
  monitors.long = c(c(61:66), c(80:85)) ,
  stimulus.start = data.frame(
    day = 4,
    hour = 10,
    minute = 04,
    second = 26
  ),
  stimulus.end = data.frame(
    day = 4,
    hour = 10,
    minute = 04,
    second = 53
  ),
  wd2 = "c1d2"
)
day2list = list(
  monitors.long = c(c(61:66), c(80:85)) ,
  stimulus.start = data.frame(
    day = 5,
    hour = 10,
    minute = 11,
    second = 46
  ),
  stimulus.end = data.frame(
    day = 5,
    hour = 10,
    minute = 12,
    second = 13
  ),
  wd2 = "c1d3"
)


day3list = list(
  monitors.long = c(c(61:66), c(80:85)) ,
  stimulus.start = data.frame(
    day = 6,
    hour = 10,
    minute = 24,
    second = 27
  ),
  stimulus.end = data.frame(
    day = 6,
    hour = 10,
    minute = 24,
    second = 54
  ),
  wd2 = "c1d4"
)
day4list = list(
  monitors.long = c(c(61:66), c(80:85)) ,
  stimulus.start = data.frame(
    day = 7,
    hour = 10,
    minute = 20,
    second = 28
  ),
  stimulus.end = data.frame(
    day = 7,
    hour = 10,
    minute = 20,
    second = 55
  ),
  wd2 = "c1d5"
)
day5list = list(
  monitors.long = c(c(61:66), c(80:85)) ,
  stimulus.start = data.frame(
    day = 8,
    hour = 09,
    minute = 56,
    second = 15
  ),
  stimulus.end = data.frame(
    day = 8,
    hour = 09,
    minute = 56,
    second = 56
  ),
  wd2 = "c1d6"
)


ref.list = list(
  day1list,
  day2list,
  day3list,
  day4list,
  day5list
)
######################3
###dayly analysis###
#####################
#this is the part that eventually I want to loop over
out.days = foreach(f = c(1:5), .combine = "rbind") %do% {
  
  # f   = 4
  ref.values = ref.list[[f]]
  wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/week2crossdata/", ref.values[[4]]  )
  setwd( wd )
  #now use monitor list to read in files
  monitor.vector = ref.values[[1]]
  #use a foreach loop to load in data perform the same data fixing to each data table
  # I found the issue! the presence of multiple days is throwing off the cmin column
  out = foreach(h = c(1:length(monitor.vector))) %do% {
    # h = 1
    Mdata <- fread(paste0("./Monitor", monitor.vector[h], ".txt")) #every well is full
    Mdata$monitor = h
    
    
    #split time
    Mdata[,c("hour","minute","second"):= tstrsplit(V3, ":")]
    Mdata[,c("day", "month", "year"):= tstrsplit(V2, "\\ ")]
    
    #how we want to lay out the data: we want to have all time seen as relative to the startle event. listed as time post startle (positive valeus) or pre startle(negavite values)
    
    #remove first 10 columns
    filtered = Mdata[,-c(1:10)]
    #find time of startle end- in minutes
    startle.end = ref.values[[3]]
    
    s.end = (as.numeric(startle.end$day) * 1440 + as.numeric(startle.end$hour) * 60) + as.numeric(startle.end$minute) + as.numeric(startle.end$second)/60
    #find cumulative time passed for each 
    filtered$cmin = (as.numeric(filtered$day) * 1440 + as.numeric(filtered$hour) * 60) + as.numeric(filtered$minute) + as.numeric(filtered$second)/60
    #subtract s.end to find time relative to startle 
    filtered$cmin = filtered$cmin -  s.end
    
    #remove unnecessary time variables
    smaller = filtered %>% 
      dplyr::select(!c(hour, minute, second, day, month, year))
    #now melt data with sex, monitor, and cmin as id variables, all other columns as variables
    meltdata = melt(smaller, id.vars = c("cmin", "monitor"), variable.name = "well.number", value.name = "activity")
    meltdata
  }
  
  bind.data = rbindlist(out)
  
  sumdata = bind.data 
  #next step is to bring in the fly genotype information, and merge taht in based on wells. have to change well info to match what is in our excel files
  #oh my gosh- one last issue, the top 6 monitors are refered to as 7-12 (in the same order)
  sumdata$well.number = gsub("V", "", sumdata$well.number) #remove V, then subtract 10
  sumdata$well.number = as.numeric(sumdata$well.number) - 10
  sumdata$wells = paste0("M.", sumdata$monitor, ".", sumdata$well.number)
  
  
  # NEED TO ADD 10 TO MATCH TRIKINETIC NOTATION
  #for well guide 4- had to load in and make changes 
  flygeno.key = fread(paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/wellguides/well.guide2." , f, ".csv"), header = T)
  #have to swap k5b-ds3 with k2b-ds3
  
  flygeno.key = flygeno.key[,-1]#first column is uneeded
  #merge in fly ids by well
  mergedata = as.data.table( merge(sumdata, flygeno.key, by = "wells"))
  colnames(mergedata)[6] = "fly.geno"
  #remove time points that are within the startle
  startle.start = ref.values[[2]]
  s.start =  (as.numeric(startle.start$day) * 1440 + as.numeric(startle.start$hour) * 60) + as.numeric(startle.start$minute) + as.numeric(startle.start$second)/60
  s.start = s.start - s.end
  #remove any time points in between 0 (s.end) and the s.start
  mergedata = mergedata[cmin < (s.start-1) | cmin > 0 ]
  
  #now find the cumulative activity for each fly for the hour prior to startle end.
  pre.data = mergedata %>% 
    group_by(fly.geno) %>% 
    filter(cmin < -2 & cmin >= -62) %>% 
    summarise(basal_act = sum(activity)) %>% 
    as.data.table(.)
  #data cleaning- we clean out any wells marked as empty- as well as any wells with basal activity of zero.
  #these are wells that either did not recieve the fly intended for it, or that fly died prior this measurement.
  mergedata = mergedata[fly.geno != "empty"]
  final = merge(mergedata, pre.data, by = "fly.geno")
  final = final[basal_act > 0 ]
  #last bit of analysis before we bind- add in day
  final$day = f
  final$fly.id = paste0(final$fly.geno, "_day", f)
  
  final
  
} #end of daily loop

wd =  paste0("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
setwd( wd )
saveRDS(out.days, paste0("week2.def.data"))
