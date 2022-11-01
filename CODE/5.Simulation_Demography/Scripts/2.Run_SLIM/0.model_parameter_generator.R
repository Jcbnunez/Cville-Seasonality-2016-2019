# Create parameter files to test slimulations
# Connor Murray 10.31.2022
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Library
library(data.table)
library(tidyverse)

# Exponential start and constant growth population
StartPop = 1500 # Should equal burn-in VCF samples
nMax = unique(c(100,250,500,750,
                         seq(from=1000, to=10000, by=1000),
                         seq(from=10000, to=50000, by=1000))) # Maximum population size

bottleneck = c(0) # Magnitude of bottleneck %
nMin = unique(c(100,250,500,750,
                         seq(from=1000, to=10000, by=1000),
                         seq(from=10000, to=50000, by=1000)))

Rep = 50 # Number of VCFs to output (from unique simulation seeds)
nSamp = 50 # Number of samples for output VCF
Gen = 50 # End of simulation

# Expanding parameters
dt <- as.data.table(expand.grid(nMax, nMin, Rep, nSamp, Gen))
setnames(dt, names(dt), c("nMax", "nMin", "Rep", "nSamp", "Gen"))

# Create variables
dt <- data.table(cbind(data.table(slurmID=c(1:c(dim(dt)[1]))), dt) %>% 
                   group_by(nMin, nMax) %>% 
                   mutate(rowID=row_number(),
                          run=paste(nMax,nMin,sep="_")))

# Reorder
setcolorder(dt, c("slurmID", "nMax", "nMin", "Rep", "nSamp", "Gen", "rowID", "run"))

# Remove lines with nMin > nMax
dt <- dt[!nMin > nMax]

# Show models being tested
ggplot(dt, 
   aes(x=as.integer(nMin), 
       y=as.integer(nMax))) +
  geom_tile(color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=-40))

# Previous run
mods <- data.table(fread("/scratch/csm6hg/bottleneck/ran.models", header = F) %>% 
                     mutate(nMax=tstrsplit(V1, "_")[[1]],
                            nMin=tstrsplit(V1, "_")[[2]]))

# Restrict to models not run
dt <- dt[!run %in% unique(mods$V1)]

# Show models being tested
ggplot() +
  geom_point(data=dt[nMin<=50000][nMax<=50000], 
         aes(x=as.integer(nMin), y=as.integer(nMax)), color="blue") +
  geom_point(data=mods[nMin<=50000][nMax<=50000], 
         aes(x=as.integer(nMin), y=as.integer(nMax)), color="red") +
  theme(axis.text.x = element_text(angle=-40))

# Remove previous runs
#dt1 <- dt[!c(nMin %in% unique(mods$nMin) & nMax %in% unique(mods$nMax))]

# Fix slurmid column
dt1 <- data.table(dt %>% mutate(slurmID=c(1:c(dim(dt)[1]))))

# Show models being tested
ggplot(dt1, 
       aes(x=as.factor(nMin), 
           y=as.factor(nMax))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle=-40))

# Output parameter file
write.csv(dt1 %>% dplyr::select(-c(run)), 
          quote=F, 
          row.names=F,
          file="/scratch/csm6hg/bottleneck/model_paramList_fin")

# Show models being tested
ggplot(dt, 
       aes(x=as.integer(nMin), 
           y=as.integer(nMax))) +
  geom_tile(color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=-40))

