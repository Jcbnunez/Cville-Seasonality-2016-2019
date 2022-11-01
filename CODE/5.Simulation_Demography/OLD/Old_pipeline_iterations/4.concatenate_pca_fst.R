# Compiles and analyze genomic data from SLiM and PLINK
# Connor Murray 11.16.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 6
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(doParallel)
library(viridis)

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- arg[1] # pathway to output
file <- arg[2] # seed for vcf
out.name <- arg[3] # the output file name

path.name="/project/berglandlab/connor/BACKUP_scratch/slim_bottleneck/data_new2"
file='21772.csv'
out.name="/scratch/csm6hg/slim_bottleneck/data/freq.out"

# VCF output folder
setwd(path.name)

# All frequency data
filenames <- list.files(pattern = file)
filenames1 <- filenames[!filenames %like% "pca"]

# Register cores
# doParallel::registerDoParallel(cores = 3)

# ReadVCF forloop - calculate diversity statistics.
out <- foreach(i = 1:length(filenames1), .combine="rbind", .errorhandling="remove") %do% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(filenames1))*100, digits = 3), 
              sep=" "))
  
  # Load allele frequency information
  dt <- data.table(read.csv(filenames1[i]),
                   file.number=i)
              
  # Finish
  return(dt)
}

# Add within years
out2 <- data.table(tot %>% 
    mutate(d.sim=abs(sim.gen2-sim.gen1)) %>% 
    mutate(dd.sim = case_when(sim.gen1 %in% c(1:16) & sim.gen2 %in% c(1:16) ~ "1.within", 
                              sim.gen1 %in% c(1:33) & sim.gen2 %in% c(17:33) | sim.gen1 %in% c(17:33) & sim.gen2 %in% c(1:33) ~ "2.overwinter", 
                              sim.gen1 %in% c(1:50) & sim.gen2 %in% c(34:50) | sim.gen1 %in% c(34:50) & sim.gen2 %in% c(1:50) ~ "3.multi")))

# Confidence interval functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

out2 %>% 
  group_by(nMax,nMin,replicate) %>% 
  summarise(length(replicate))

# Average across replicates
out.fin <- data.table(out2 %>% 
                        group_by(nMax, nMin, dd.sim) %>% 
                        summarize(mean.FST = mean(FST, na.rm = T),
                                  sd.FST = sd(FST, na.rm = T),
                                  n.FST = length(unique(replicate))) %>% 
                        mutate(se.FST = sd.FST/sqrt(n.FST)) %>% 
                        mutate(lci.FST = lower_ci(mean.FST, se.FST, n.FST, conf_level = 0.95),
                               uci.FST = upper_ci(mean.FST, se.FST, n.FST, conf_level = 0.95)))

# Output aggregated FST
write.csv(out.fin, file=paste("../data.fst.aggregate.", 
                              unique(out.fin$nMax), 
                              unique(out.fin$nMin)), 
          quote = F, row.names = F)

# ReadVCF forloop - calculate diversity statistics.
out2 <- foreach(i = 15173:length(filenames1), .combine="rbind", .errorhandling="remove") %do% {
  
  # Progress message
  print(paste("Progress:", 
              round((i/length(filenames1))*100, digits = 3), 
              sep=" "))
  
  # Load allele frequency information
  dt <- data.table(read.csv(filenames1[i]),
                   file.number=i)
  
  # Finish
  return(dt)
}

# write.csv(out2, file="../data.fst.100reps.more2.csv", quote = F, row.names = F)

# Read in FST info
tot <- data.table(rbind(read.csv("../data.fst.100reps.more1.csv", 
                                 header = T),
                        read.csv("../data.fst.100reps.more2.csv", 
                                 header = T)))

# Add within years
out2 <- data.table(tot %>% 
  mutate(d.sim=abs(sim.gen2-sim.gen1)) %>% 
  mutate(dd.sim = case_when(sim.gen1 %in% c(1:16) & sim.gen2 %in% c(1:16) ~ "1.within", 
                            sim.gen1 %in% c(1:33) & sim.gen2 %in% c(17:33) | sim.gen1 %in% c(17:33) & sim.gen2 %in% c(1:33) ~ "2.overwinter", 
                            sim.gen1 %in% c(1:50) & sim.gen2 %in% c(34:50) | sim.gen1 %in% c(34:50) & sim.gen2 %in% c(1:50) ~ "3.multi")))

# Some failed during pairwise FST calculation
tottest <- data.table(tot %>% 
                        group_by(nMax, nMin, replicate, seed) %>% 
                        summarise(replicate.num=length(replicate)))

# Remake
redo <- tottest[!replicate.num %in% 1225] %>% 
  group_by(nMax, nMin) %>% 
  summarise(unique(nMax),
            unique(nMin))

write.csv(redo, "/scratch/csm6hg/slim_bottleneck/redo.sims.csv", quote = F, row.names = F)

# Confidence interval functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

out2 %>% 
  group_by(nMax,nMin,replicate) %>% 
  summarise(length(replicate))

# Average across replicates
out.fin <- data.table(out2 %>% 
            group_by(nMax, nMin, dd.sim) %>% 
            summarize(mean.FST = mean(FST, na.rm = T),
                      sd.FST = sd(FST, na.rm = T),
                      n.FST = length(unique(replicate))) %>% 
            mutate(se.FST = sd.FST/sqrt(n.FST)) %>% 
            mutate(lci.FST = lower_ci(mean.FST, se.FST, n.FST, conf_level = 0.95),
                   uci.FST = upper_ci(mean.FST, se.FST, n.FST, conf_level = 0.95)))

# Output aggregated FST
write.csv(out.fin, file="../data.fst.aggregate.100reps.more.csv", quote = F, row.names = F)

# Add comparisons
out.fin <- data.table(out.fin %>% mutate(comp=paste(nMax,nMin,sep="_")))

k <- data.table(comparison=unique(out.fin$comp),
                num=unique(out.fin$n.FST)) %>% 
      mutate(num.left=100-num)







pdf("../FST.pairwise.ribbon.100reps.1k-100k.more.pdf", width = 13, height = 8)

# Plot within vs yearly FST
ggplot(out.fin,
         aes(x=as.numeric(nMin),
             y=as.numeric(mean.FST),
             ymin=as.numeric(lci.FST),
             ymax=as.numeric(uci.FST),
             color=dd.sim)) +
    geom_ribbon(alpha=0.5, 
                aes(fill=dd.sim)) +
    geom_line() +
    geom_hline(yintercept = 0.0, linetype=2, size=1.1) +
    scale_x_log10() +
    facet_wrap(~nMax) +
    labs(x="Minimum", 
         y="Mean FST", 
         fill="Year") +
    theme_classic() +
    guides(color=FALSE) +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=12),
          legend.title = element_text(face="bold", size=15),
          legend.background = element_blank(),
          strip.text =element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=18, 
                                     angle = -30, vjust = 0.6),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))

dev.off()

pdf("../FST.pairwise.raster.100reps.1k-100k.more.pdf", width = 13, height = 8)

# Plot within vs yearly FST Raster plot
ggplot(out.fin,
       aes(x=as.factor(nMin),
           y=as.factor(nMax),
           fill=mean.FST)) +
  geom_tile() +
  facet_wrap(~dd.sim, scales="free") +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="Mean FST") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18, 
                                   angle = -30, vjust = 0.6),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()

# All metapopulation data
filenames2 <- list.files(pattern = ".txt$")

# ReadVCF forloop - calculate diversity statistics.
pop <- foreach(i = 1:length(filenames2), .combine="rbind", .errorhandling="remove") %dopar% {
  
  # Progress message
  print(i)
  
  # Load allele frequency information
  dt <- data.table(fread(filenames2[i]),
                   file.number=i)
  
  # Finish
  return(dt)
}

# Column names
colnames(pop) <- c("Gen", "Fin.gen", "n", "nMax", 
                   "nMin", "seed", "slurmid", "genMax", 
                   "file.number")

# Add bottleneck
pop2 <- data.table(pop %>% 
                     group_by(seed) %>% 
                     mutate(bottle=((nMax-nMin)/nMax)*100))

pdf("../popultion.model.100reps.1k-100k.pdf", width = 13, height = 8)

# Plot within vs yearly FST
ggplot(pop2[Gen>1],
       aes(x=Gen,
           y=n,
           group=nMin,
           color=nMin)) +
  geom_line(size=1.1) +
  facet_wrap(~nMax, scales="free_x") +
  labs(x="Generation", 
       y="Population size", 
       color="Minimum population size",
       title="") +
  theme_classic() +
  #scale_y_log10() +
  scale_color_viridis(option = "magma", begin = 1, end = 0) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=13),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

dev.off()

# All frequency data
filenames3 <- list.files(pattern = "pca")

# ReadVCF forloop - calculate diversity statistics.
out.pca <- foreach(i = 1:length(filenames3), .combine="rbind", .errorhandling="remove") %dopar% {

  # Progress message
  print(paste("Progress:", 
              round((i/length(filenames3))*100, digits = 3), 
              sep=" "))
  
  # Load allele frequency information
  dt <- data.table(read.csv(filenames3[i]),
                   r2pc1=summary(lm(data = read.csv(filenames3[i]), Dim.1 ~ Year))$r.squared,
                   r2pc2=summary(lm(data = read.csv(filenames3[i]), Dim.2 ~ Year))$r.squared,
                   r2pc3=summary(lm(data = read.csv(filenames3[i]), Dim.3 ~ Year))$r.squared,
                   file.number=i,
                   nMin=tstrsplit(filenames3[i], ".", fixed=T)[[5]],
                   nMax=tstrsplit(filenames3[i], ".", fixed=T)[[4]],
                   replicate=tstrsplit(filenames3[i], ".", fixed=T)[[6]],
                   seed=tstrsplit(filenames3[i], ".", fixed=T)[[7]])
  
  # Finish
  return(dt)
}

# Output pca
#write.csv(out.pca, file="../data.pca.100reps.more.csv", quote = F, row.names = F)

# Average across replicates
out.fin <- data.table(out.pca %>% 
                        group_by(nMax, nMin, Year) %>% 
                        summarize(mean.pca1 = mean(Dim.1, na.rm = T),
                                  mean.pca2 = mean(Dim.2, na.rm = T),
                                  mean.pca3 = mean(Dim.3, na.rm = T),
                                  mean.r2pc1 = mean(r2pc1, na.rm = T),
                                  mean.r2pc2 = mean(r2pc2, na.rm = T),
                                  mean.r2pc3 = mean(r2pc3, na.rm = T),
                                  sd.pca1 = sd(Dim.1, na.rm = T),
                                  sd.pca2 = sd(Dim.2, na.rm = T),
                                  sd.pca3 = sd(Dim.3, na.rm = T),
                                  sd.r2pc1 = sd(r2pc1, na.rm = T),
                                  sd.r2pc2 = sd(r2pc2, na.rm = T),
                                  sd.r2pc3 = sd(r2pc3, na.rm = T),
                                  n.pca = length(unique(replicate))) %>% 
                        mutate(se.pca1 = sd.pca1/sqrt(n.pca),
                               se.pca2 = sd.pca2/sqrt(n.pca),
                               se.pca3 = sd.pca3/sqrt(n.pca),
                               se.r2pc1 = sd.r2pc1/sqrt(n.pca),
                               se.r2pc2 = sd.r2pc2/sqrt(n.pca),
                               se.r2pc3 = sd.r2pc3/sqrt(n.pca)) %>% 
                        mutate(lci.pca1 = lower_ci(mean.pca1, se.pca1, n.pca, conf_level = 0.95),
                               uci.pca1 = upper_ci(mean.pca1, se.pca1, n.pca, conf_level = 0.95),
                               lci.pca2 = lower_ci(mean.pca2, se.pca2, n.pca, conf_level = 0.95),
                               uci.pca2 = upper_ci(mean.pca2, se.pca2, n.pca, conf_level = 0.95),
                               lci.pca3 = lower_ci(mean.pca3, se.pca3, n.pca, conf_level = 0.95),
                               uci.pca3 = upper_ci(mean.pca3, se.pca3, n.pca, conf_level = 0.95),
                               lci.r2pc1 = lower_ci(mean.r2pc1, se.r2pc1, n.pca, conf_level = 0.95),
                               uci.r2pc1 = upper_ci(mean.r2pc1, se.r2pc1, n.pca, conf_level = 0.95),
                               lci.r2pc2 = lower_ci(mean.r2pc2, se.r2pc2, n.pca, conf_level = 0.95),
                               uci.r2pc2 = upper_ci(mean.r2pc2, se.r2pc2, n.pca, conf_level = 0.95),
                               lci.r2pc3 = lower_ci(mean.r2pc3, se.r2pc3, n.pca, conf_level = 0.95),
                               uci.r2pc3 = upper_ci(mean.r2pc3, se.r2pc3, n.pca, conf_level = 0.95)))

# Output aggregated pca
#write.csv(out.fin, file="../data.pca.aggregate.100reps.more.csv", quote = F, row.names = F)

pdf("../pca.raster.100reps.1k-100k.more.pdf", width = 13, height = 8)

# Plot within vs yearly FST Raster plot
ggplot(out.pca,
       aes(x=reorder(as.factor(nMin), as.numeric(nMin)),
           y=reorder(as.factor(nMax), as.numeric(nMax)),
           fill=as.numeric(r2pc1))) +
  geom_raster() +
  #facet_wrap(~Year) +
  labs(x="Minimum population size", 
       y="Maximum population size", 
       fill="R2 PC1 ~ Year") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=18, 
                                   angle=-30, vjust=0.6),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()
