# Create parameter files to test slimulations
# Connor Murray 10.19.2021
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Library
library(data.table)
library(tidyverse)

# Exponential start and constant growth population
StartPop = 1500 # Should equal burn-in VCF samples
nMax = unique(c(50,100,500,1000,2500,7500,10000,25000,50000,75000,100000,
         seq(from=1000, to=7000, by=500),
         seq(from=10000, to=74000, by=1000))) # Maximum population size
bottleneck = c(0) # Magnitude of bottleneck %
nMin = unique(c(nMax, 50,100,500,1000,2500,7500,10000,25000,50000,75000,100000, 
         seq(from=1000, to=7000, by=500)))
Rep = 100 # Number of VCFs to output (from unique simulation seeds)
nSamp = 50 # Number of samples for output VCF
Gen = 50 # End of simulation

# Expanding parameters
dt <- as.data.table(expand.grid(nMax, nMin, Rep, nSamp, Gen))
setnames(dt, names(dt), c("nMax", "nMin", "Rep", "nSamp", "Gen"))

# Create variables
dt <- data.table(cbind(data.table(slurmID=c(1:c(dim(dt)[1]))), dt) %>% 
                   group_by(nMin, nMax) %>% 
                   mutate(rowID=row_number()))

# Reorder
setcolorder(dt, c("slurmID", "nMax", "nMin", "Rep", "nSamp", "Gen", "rowID"))

# Remove lines with nMin > nMax
dt <- dt[!nMin > nMax]

# Show models being tested
ggplot(dt, 
   aes(x=as.factor(nMin), 
       y=as.factor(nMax))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle=-40))

# Fix slurmid column
dt <- data.table(dt %>% mutate(slurmID=c(1:c(dim(dt)[1]))))

# Output parameter file
write.csv(dt, quote=F, row.names=F,
         file="/scratch/csm6hg/slim_bottleneck/model_paramList6")
