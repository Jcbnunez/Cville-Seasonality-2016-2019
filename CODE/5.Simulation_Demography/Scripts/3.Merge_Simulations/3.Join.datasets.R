# Find and merge output data from SLIM
# Connor S. Murray
# 10.31.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3;module load gdal geos proj; R

library(foreach)
library(data.table)
library(tidyverse)
library(foreach)

####
root <- "/scratch/csm6hg/bottleneck/Parsed_sims"
setwd(root)
files <- system( paste("ls ", root, " | grep '.Rdata'  "), intern = T )
joint.dat <- foreach(i=1:length(files), .combine = "rbind", .errorhandling = "pass") %do% {
  message(paste(i, "/",length(files)))
  tmp <- get(load( paste( files[i], sep = "" ) ))
}

# Save output
saveRDS("/project/berglandlab/connor/bottleneck/parsed.data.Rdata", object = joint.dat)

pp <- data.table(joint.dat %>% group_by(model) %>% 
                   summarize(num.reps=table(model)))
ran.models <- unique(pp[num.reps.N>=50]$model)

# Output
write.table(x = ran.models, file = "/scratch/csm6hg/bottleneck/ran.models", 
            quote = F, row.names = F, col.names = F)

# Restrict to analyses 
restrict <- data.table(fread("../model_paramList_fin4") %>% 
                         mutate(model=paste(nMax, nMin, sep="_")))

# Restrict simulations to runs
joint.dat <- joint.dat[nMax<=50000][nMin<=50000][nMax %in% 
              c(250,500,750, seq(from=1000, to=50000, by=1000))][nMin %in% 
              c(250,500,750, seq(from=1000, to=50000, by=1000))]

# Show sims
ggplot(joint.dat, aes(x=as.numeric(nMin), y=as.numeric(nMax))) + geom_point()
