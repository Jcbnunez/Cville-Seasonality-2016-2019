# Compiles and analyze genomic data from SLiM output
# Connor Murray 7.20.2022
# ijob -A berglandlab --mem=50G -p largemem -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(doParallel)
library(viridis)

# VCF output folder
setwd("/project/berglandlab/connor/bottleneck/data3")

# All frequency data
filenames1 <- list.files(pattern = ".csv$")
file.pca <- filenames1[filenames1 %like% "pca"]
file.afvar <- filenames1[filenames1 %like% "afvar"]
file.euc <- filenames1[filenames1 %like% "euc"]
file.fst <- filenames1[filenames1 %like% "fst"]

# Register cores
#doParallel::registerDoParallel(cores = 10)

# Extract metadata from filenames
meta <- data.table(data.table(seed=tstrsplit(file.pca, ".", fixed=T)[[8]],
                              nMax = tstrsplit(file.pca, ".", fixed=T)[[5]],
                              nMin = tstrsplit(file.pca, ".", fixed=T)[[6]],
                              replicate = tstrsplit(file.pca, ".", fixed=T)[[7]],
                              file=file.pca) %>% 
                     group_by(nMax, nMin) %>% 
                     mutate(run_number=row_number(),
                            run=paste(nMax, nMin,sep="_")))

# Extract metadata from filenames
meta.fst <- data.table(data.table(seed=tstrsplit(file.pca, ".", fixed=T)[[8]],
                              nMax = tstrsplit(file.pca, ".", fixed=T)[[5]],
                              nMin = tstrsplit(file.pca, ".", fixed=T)[[6]],
                              replicate = tstrsplit(file.pca, ".", fixed=T)[[7]],
                              file=file.pca) %>% 
                     group_by(nMax, nMin) %>% 
                     mutate(run_number=row_number(),
                            run=paste(nMax, nMin,sep="_")))


# Distinct models
runs <- unique(meta$run)

# Restrict to models we want
mods <- data.table(read.csv("/scratch/csm6hg/bottleneck/model_paramList_fin3") %>% 
                     mutate(run=paste(nMax,nMin, sep="_")) %>% 
                     distinct(nMax, nMin, run))

# Restrict to models
meta <- data.table(meta[!run %in% mods$run] %>% 
  group_by(run) %>% 
  mutate(num_runs=table(run)))

meta <- meta[num_runs.N==100]

# Load Fst information
dt2 <- data.table(do.call(rbind, lapply(file.fst, fread)))

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

# Plot yearly change in FST
fst.plot <- {dt2[!nMax<5000] %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("mean.FST_Within",
                        "mean.FST_Overwinter",
                        "mean.FST_Multi")) %>%
    mutate(Year=case_when(name=="mean.FST_Within" ~ 1,
                          name=="mean.FST_Overwinter" ~ 2,
                          name=="mean.FST_Multi" ~ 3)) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    ggplot(., 
           aes(x=as.factor(Year), 
               y=value, 
               group=interaction(replicate, model),
               color=log10(nMin))) + 
    geom_point(size=2.2, alpha=0.6) +
    geom_line(size=1, alpha=0.6) +
    facet_wrap(~bot.class, nrow=1) +
    theme_classic() +
    labs(x="Number of Winters (sim. years)",
         y="Fst",
         color="log10(Minimum pop.)") +
    scale_color_viridis(begin = 1, end=0) +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

# ReadVCF forloop - PCA and year correlation
out.pca <- foreach(i = 1:length(runs), .combine="rbind", .errorhandling="remove") %do% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(runs))*100, digits = 3), 
              sep=" "))
  
  # All replicate files per run
  filenames4 <- meta[run==runs[i]]$file
  
  # Load allele frequency information
  dt.i <- data.table(do.call(rbind, lapply(filenames4, fread)),
                   run.number=i,
                   run=runs[i]) %>% 
    mutate(nMax=as.integer(tstrsplit(Sample, "_")[[7]]),
           nMin=as.integer(tstrsplit(Sample, "_")[[8]]),
           replicate=as.integer(tstrsplit(Sample, "_")[[9]]),
           seed=tstrsplit(Sample, "_")[[10]])
  
  # dt.i[replicate==1] %>% ggplot(., aes(x=Dim.1, y=Dim.2, color=Year)) + geom_point()
  # dt.i[replicate==1] %>% ggplot(., aes(x=LD1, y=LD2, color=Year)) + geom_point()
  
  # Aggregate simulation mean FST and PCA 
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
  dt.i.fin <- data.table(dt.r2 %>% 
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
  
  # Finish
  return(dt.i.fin)
}

# ReadVCF forloop - calculate diversity statistics.
out.pca1 <- foreach(i = 1:length(runs), .combine="rbind", .errorhandling="remove") %do% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(runs))*100, digits = 3), 
              sep=" "))
  
  # All replicate files per run
  filenames4 <- meta[run==runs[i]]$file
  
  # Load allele frequency information
  dt.i <- data.table(do.call(rbind, lapply(filenames4, fread)),
                     run.number=i,
                     run=runs[i]) %>% 
    mutate(nMax=as.numeric(tstrsplit(Sample, "_")[[7]]),
           nMin=as.numeric(tstrsplit(Sample, "_")[[8]]),
           replicate=tstrsplit(Sample, "_")[[9]],
           seed=tstrsplit(Sample, "_")[[10]]) %>% 
    mutate(bottle=(nMax-nMin)/nMax) %>% 
    mutate(bot.class=case_when(bottle >=0.9 ~ "Perc90",
                               bottle >=0.8 & bottle < 0.9 ~ "Perc80-90",
                               bottle >= 0.5 & bottle < 0.8 ~ "Perc50-80",
                               bottle < 0.5 ~ "Less50"))
  
  # dt.i[replicate==1] %>% ggplot(., aes(x=Dim.1, y=Dim.2, color=Year)) + geom_point()
  # dt.i[replicate==1] %>% ggplot(., aes(x=LD1, y=LD2, color=Year)) + geom_point()

}

# Load allele frequency variance information
dt.var <- data.table(do.call(rbind, lapply(file.afvar, fread)))

# Average across replicates
dt2.var <- data.table(dt.var %>% 
                    group_by(nMax, nMin, model, replicate) %>% 
                    summarize(mean.var = mean(Mean.var, na.rm = T),
                              median.var = mean(Median.var),
                              IQR05.var = mean(IQR05.var),
                              IQR95.var = mean(IQR95.var)))

# Load euclidean distance information
dt.euc <- data.table(do.call(rbind, lapply(file.euc, fread)))

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
                merge(out.pca, 
                         by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                merge(dt2.var, 
                         by=c('model', 'nMax', 'nMin', 'replicate')) %>% 
                merge(dt2.euc, 
                         by=c('model', 'nMax', 'nMin', 'replicate')))

# Output aggregated file
write.csv(mean.sim, file="../data.sim.100reps.csv", quote = F, row.names = F)
