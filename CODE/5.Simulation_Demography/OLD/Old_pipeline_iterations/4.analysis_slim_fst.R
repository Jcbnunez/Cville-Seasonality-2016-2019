# Compiles and analyze genomic data from SLiM and PLINK
# Connor Murray 10.13.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 1
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(viridis)

# VCF output folder
setwd("/scratch/csm6hg/slim_bottleneck/data")

# All frequency data
filenames1 <- list.files(pattern = ".csv$")

# ReadVCF forloop - calculate diversity statistics.
out <- foreach(i = 1:length(filenames1), .combine="rbind", .errorhandling="remove") %do% {
  
  # Progress message
  print(paste("Progress:", 
        round((i/length(filenmes1))*100, digits = 3), 
        sep=" "))
  
  # Load allele frequency information
  dt <- data.table(read.csv(filenames1[i]),
                   file.number=i)
              
  # Finish
  return(dt)
}

# Add within years
out2 <- out %>% 
  mutate(d.sim=abs(sim.gen2-sim.gen1)) %>% 
  mutate(dd.sim = case_when(d.sim <=14 & d.sim > 3 ~ "1.within",
                            d.sim >=16 & d.sim <=33 ~ "2.overwinter",
                            d.sim >= 34 ~ "3.multi")) 

# Average across replicates
out.fin <- data.table(out2 %>% 
             group_by(nMax, nMin, bottle, dd.sim) %>% 
            summarize(mean.FST=mean(FST, na.rm = T),
                      lci.FST=mean(FST, na.rm = T) - 1.96*(sd(FST)/sqrt(100)),
                      uci.FST=mean(FST, na.rm = T) + 1.96*(sd(FST)/sqrt(100))))

pdf("../FST.pairwise.pdf", width = 13, height = 8)

# Plot within vs yearly FST
ggplot(na.omit(out.fin),
         aes(x=nMin,
             y=mean.FST,
             ymin=lci.FST,
             ymax=uci.FST,
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
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))
  
dev.off()

# All metapopulation data
filenames2 <- list.files(pattern = ".txt$")

# ReadVCF forloop - calculate diversity statistics.
pop <- foreach(i = 1:length(filenames2), .combine="rbind", .errorhandling="remove") %do% {
  
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

pdf("../popultion.model.pdf", width = 13, height = 8)

# Plot within vs yearly FST
ggplot(pop2,
       aes(x=Gen,
           y=n,
           group=bottle,
           color=bottle)) +
  geom_line(size=1.1) +
  facet_wrap(~nMax, scales="free_x") +
  labs(x="Generation", 
       y="Population size", 
       color="Bottleneck size (%)",
       title="") +
  theme_classic() +
  #scale_y_log10() +
  scale_color_viridis(option = "magma", begin = 1, end = 0) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold.italic", size=10),
        legend.title = element_text(face="bold", size=13),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=15),
        axis.text.y = element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

dev.off()