# Create parameter files to test slimulations
# Connor Murray 10.19.2021
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Library
library(data.table)
library(tidyverse)

# Exponential start and constant growth population
StartPop = 1500 # Should equal burn-in VCF samples
nMax = seq(from=10000, to=74000, by=1000) # Maximum population size
#bottleneck = c(0, 0.25, 0.5, 0.8, seq(0.9, 0.99, by=0.01)) # Magnitude of bottleneck (%)
bottleneck = c(0) # Magnitude of bottleneck %
nMin = seq(from=1000, to=7000, by=100)
Rep = 100 # Number of VCFs to output (from unique simulation seeds)
nSamp = 50 # Number of samples for output VCF
Gen = 50 # End of simulation

# Expanding parameters
dt <- as.data.table(expand.grid(nMax, nMin, Rep, nSamp, Gen))
setnames(dt, names(dt), c("nMax", "nMin", "Rep", "nSamp", "Gen"))

# Create variables
dt <- data.table(cbind(data.table(slurmID=c(1:c(dim(dt)[1]*Rep))), dt) %>% 
                   mutate(simID=round(abs(rnorm(1, 1e8, n=c(1:dim(dt)[1]*Rep))), digits=0)) %>% 
                   group_by(nMin, nMax) %>% 
                   mutate(rowID=row_number()))

# Reorder
setcolorder(dt, c("slurmID", "nMax", "nMin", "Rep", "nSamp", "Gen", "simID"))

# Remove lines with nMin > nMax
dt <- dt[!nMin > nMax]

# Fix slurmid column
dt <- data.table(dt %>% mutate(slurmID=c(1:c(dim(dt)[1]))))

# Exclude already ran samples
dt2 <- data.table(read.csv("/scratch/csm6hg/slim_bottleneck/model_paramList3"))

dt[nMax %in% dt2$nMax & 
   nMin %in% dt2$nMin]

# Output parameter file
write.csv(dt, quote=F, row.names=F,
         file="/scratch/csm6hg/slim_bottleneck/model_paramList4")
