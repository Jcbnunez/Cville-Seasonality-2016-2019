# Compiles and analyze genomic data from SLiM output
# Connor Murray 7.27.2022
# ijob -A berglandlab --mem=50G -p largemem -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(doParallel)
library(viridis)

# Register cores
#doParallel::registerDoParallel(cores = 10)

# VCF output folder
setwd("/scratch/csm6hg/bottleneck/data.condense")

# All frequency data
filenames1 <- list.files(pattern = ".csv$")

# Extract metadata from filenames
meta <- data.table(data.table(nMax = rev(tstrsplit(filenames1, ".", fixed=T))[[3]],
                              nMin = rev(tstrsplit(filenames1, ".", fixed=T))[[2]],
                              statistic = rev(tstrsplit(filenames1, ".", fixed=T))[[4]],
                              file=filenames1) %>% 
                     group_by(nMax, nMin) %>% 
                     mutate(run_number=row_number(),
                            run=paste(nMax, nMin,sep="_")))

# Restrict to models we want
mods <- data.table(read.csv("/scratch/csm6hg/bottleneck/model_paramList_fin3") %>% 
                     mutate(run=paste(nMax, nMin, sep="_")) %>% 
                     distinct(nMax, nMin, run))

# Restrict to models
meta <- data.table(meta[run %in% mods$run] %>% 
  group_by(run) %>% 
  mutate(num_runs = table(run)))

# Distinct models
runs <- unique(meta$run)

# Other files
file.pca <- filenames1[filenames1 %like% "pca"]
file.afvar <- filenames1[filenames1 %like% "afvar"]
file.euc <- filenames1[filenames1 %like% "euc"]
file.fst <- filenames1[filenames1 %like% "fst"]

# Load PCA / LD information
dt.i <- data.table(do.call(rbind, lapply(file.pca, fread))) 
colnames(dt.i) <- c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5",
                    "pc_perc.comp1", "pc_perc.comp2", "pc_perc.comp3", 
                    "pc_perc.comp4", "pc_perc.comp5",
                    "Sample", "Sample.gen", "Year", "samp.pop",
                    "LD1", "LD2")

dt.i <- data.table(dt.i %>% 
     mutate(nMax=as.integer(tstrsplit(Sample, "_")[[7]]),
           nMin=as.integer(tstrsplit(Sample, "_")[[8]]),
           replicate=as.integer(tstrsplit(Sample, "_")[[9]]),
           seed=tstrsplit(Sample, "_")[[10]]))
  
# dt.i[replicate==1] %>% ggplot(., aes(x=Dim.1, y=Dim.2, color=Year)) + geom_point()
# dt.i[replicate==1] %>% ggplot(., aes(x=LD1, y=LD2, color=Year)) + geom_point()
  
# Aggregate simulation mean PCA 
dt.r2 <- data.table(dt.i %>% 
          mutate(dd.sim = case_when(Year %in% "Year1" ~ "1.within",
                                    Year %in% "Year2" ~ "2.overwinter",
                                     Year %in% "Year3" ~ "3.multi")) %>% 
          filter(!dd.sim == "3.multi") %>%
          group_by(replicate) %>% 
          mutate(r2dapc1=summary(lm(LD1 ~ Year))$r.squared,
                 r2dapc2=summary(lm(LD2 ~ Year))$r.squared,
                 r2pc1=summary(lm(Dim.1 ~ Year))$r.squared,
                 r2pc2=summary(lm(Dim.2 ~ Year))$r.squared,
                 r2pc3=summary(lm(Dim.3 ~ Year))$r.squared)) 
  
# Aggregate further
dt.i.pca <- data.table(dt.r2 %>%
            group_by(nMax, nMin, dd.sim, replicate) %>% 
            summarize(mean.dapc1 = mean(LD1, na.rm = TRUE),
                      mean.dapc2 = mean(LD2, na.rm = TRUE), 
                      median.dapc1 = median(LD1, na.rm = TRUE),
                      median.dapc2 = median(LD2, na.rm = TRUE), 
                      r2dapc1 = mean(r2dapc1, na.rm = TRUE),
                      r2dapc2 = mean(r2dapc2, na.rm = TRUE),
                      LD1_corr.time.est = median(r2dapc1, na.rm = TRUE),
                      LD2_corr.time.est = median(r2dapc2, na.rm = TRUE),
                      mean.pc1 = mean(Dim.1, na.rm = TRUE),
                      mean.pc2 = mean(Dim.2, na.rm = TRUE), 
                      mean.pc3 = mean(Dim.3, na.rm = TRUE),
                      r2pc1 = mean(r2pc1, na.rm = TRUE),
                      r2pc2 = mean(r2pc2, na.rm = TRUE),
                      r2pc3 = mean(r2pc3, na.rm = TRUE),
                      perc.comp1.mean = mean(pc_perc.comp1, na.rm = T),
                      perc.comp2.mean = mean(pc_perc.comp2, na.rm = T),
                      perc.comp3.mean = mean(pc_perc.comp3, na.rm = T),
                      perc.comp1.med = median(pc_perc.comp1, na.rm = T),
                      perc.comp2.med = median(pc_perc.comp2, na.rm = T),
                      perc.comp3.med = median(pc_perc.comp3, na.rm = T),
                      Dim.1_corr.time.est = median(r2pc1, na.rm = TRUE),
                      Dim.2_corr.time.est = median(r2pc2, na.rm = TRUE),
                      Dim.3_corr.time.est = median(r2pc3, na.rm = TRUE)) %>% 
            mutate(model = paste(nMax, nMin, sep="_"),
                   dd.sim = dd.sim) %>% 
            pivot_wider(names_from = dd.sim,
                        values_from=c(mean.dapc1:r2pc3)))

# Load Fst information
dt2 <- data.table(do.call(rbind, lapply(file.fst, fread)))
colnames(dt2) <- c("samp1","samp2","FST","iteration","nMax","nMin",
                   "replicate","seed","YearS1","sim.gen1","YearS2","sim.gen2","Year")

# Average across replicates
dt2 <- data.table(dt2 %>% 
                    group_by(nMax, nMin, Year, replicate) %>% 
                    summarize(mean.FST = mean(FST, na.rm = T),
                              median.FST = quantile(FST, probs = 0.5),
                              IQR05.FST = quantile(FST, probs = 0.05),
                              IQR95.FST = quantile(FST, probs = 0.95)) %>% 
                    mutate(model=paste(nMax, nMin, sep="_")) %>% 
                    pivot_wider(names_from = Year,
                                values_from=c(mean.FST:IQR95.FST)))

# Load allele frequency variance information
dt.var <- data.table(do.call(rbind, lapply(file.afvar, fread)))
colnames(dt.var) <- c("Mean.var","Median.var","IQR95.var",
                      "IQR05.var","replicate","nMax","nMin","model")

# Average across replicates
dt2.var <- data.table(dt.var %>% 
                    group_by(nMax, nMin, model, replicate) %>% 
                    summarize(mean.var = mean(Mean.var, na.rm = T),
                              median.var = mean(Median.var),
                              IQR05.var = mean(IQR05.var),
                              IQR95.var = mean(IQR95.var)))

# Load euclidean distance information
dt.euc <- data.table(do.call(rbind, lapply(file.euc, fread)))
colnames(dt.euc) <- c("year.diff","Mean.euc","Median.euc",
                      "IQR05.euc","IQR95.euc","data",
                      "replicate","nMax","nMin","model")

# Average across replicates
dt2.euc <- data.table(dt.euc %>% 
                  group_by(nMax, nMin, model, data, year.diff, replicate) %>% 
                  summarize(mean.euc = mean(Mean.euc, na.rm = T),
                            median.euc = mean(Median.euc),
                            IQR05.euc = mean(IQR05.euc),
                            IQR95.euc = mean(IQR95.euc)) %>% 
                  pivot_wider(names_from = c(data, year.diff),
                              values_from=c(mean.euc:IQR95.euc)))

# Aggregate all statistics
mean.sim <- data.table(dt2 %>% 
                merge(dt.i.pca, 
                         by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                merge(dt2.var, 
                         by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                merge(dt2.euc, 
                         by=c('model', 'nMax', 'nMin', 'replicate')))

# Output aggregated file
write.csv(mean.sim, file="../data.sim.100reps.proto.csv", quote = F, row.names = F)
