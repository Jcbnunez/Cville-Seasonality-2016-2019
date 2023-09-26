# ABC analysis and plotting of Drosophila simulations
# Connor S. Murray
# 9.14.2023
# ijob -c 4 --mem=20G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3;module load gdal geos proj; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)
library(foreach)
library(viridis)
library(abc)
library(abc.data)
library(cowplot)

# Working directory
setwd("C:/Users/Conno/Desktop/gitmaster/drosophila_bottleneck_simulations/Data")

# Simulation data
joint.dat <- readRDS("parsed.data.Rdata")

# Restrict simulations to runs below 50K
joint.dat <- joint.dat[nMax<=50000][nMin<=50000][nMax %in% 
              c(250,500,750, seq(from=1000, to=50000, by=1000))][nMin %in% 
              c(250,500,750, seq(from=1000, to=50000, by=1000))]

# Load empirical data: obs_values.final
load("obs_values.final.Rdata")

# Change names of empirical data
empir <- data.table(obs_values.final %>% 
                      pivot_wider(names_from="Obs_var", values_from="Median_Val"))
empir <- empir[,c(1:5,11,19,25,26,32:34)]
colnames(empir) <- as.character(tstrsplit(colnames(empir),"=")[[2]])

# work with sim data
mean.sim <- joint.dat
mean.sim <- mean.sim[,c(1:4,6,15:17,30,32,49:51,61:63)][nMax>=250]

# Standardize euclidean distances
mean.sim <- data.table(mean.sim %>% 
  mutate(median.euc_Euc2D_0=(median.euc_Euc2D_0/100),
         median.euc_Euc2D_1=(median.euc_Euc2D_1/100),
         median.euc_Euc2D_2=(median.euc_Euc2D_2/100)))

colnames(mean.sim)[c(5,9:16)] <- c("50%_AF_var", "LD1_corr.time.est", "LD2_corr.time.est", 
                                   "3.Multi-Year_FST", "2.Overwinter_FST", "1.within_FST",
                                   "0_Median.standarized.PCAEucDist",
                                   "1_Median.standarized.PCAEucDist", 
                                   "2_Median.standarized.PCAEucDist")

# Model information
models <- data.table(nMax=as.numeric(tstrsplit(mean.sim$model, "_")[[1]]),
                     nMin=as.numeric(tstrsplit(mean.sim$model, "_")[[2]])) %>% 
  mutate(bottle=((nMax-nMin)/nMax)*100)

# Statistics to use
sstats <-  c(colnames(empir))
sstats <- sstats[-c(2,3,6,7,8,10)]

# Empirical stats
empir <- empir %>% dplyr::select(sstats)

# ABC stats
sim.stat <- mean.sim %>% dplyr::select(sstats) %>% 
  relocate(colnames(empir))

# Thresholds to test in ABC
thresh <- seq(from=0.01, to = 0.15, by = 0.01)

# Go through several thresholds
tmp.out <- foreach(i=1:length(thresh), .combine = "rbind") %do% {
  # i=5
  
  # ABC using rejection methods
  tmp <- abc(target = empir, 
             param = models %>% dplyr::select(-c(bottle)),
             sumstat = sim.stat,
             tol=thresh[i], 
             transf=c("none","log"),
             method="loclinear")
  
  # Posterior probabilities - daignostic plot
  # plot(tmp, models %>% dplyr::select(-c(bottle)), file="diagnostic.plots.abc")
  
  # ESS of posterior distribution
  ess.post <- apply(tmp$adj.values, 2, coda::effectiveSize)
  
  # Compile data
  tmp <- data.table(summary(tmp), 
                    empir,
                    threshold=thresh[i],
                    ess.post)
  
  # Progress message
  print(paste("Threshold:", thresh[i], sep=" "))
  
  # Finish
  return(tmp)
}

# Fix data
colnames(tmp.out)[1:3] <- c("data", "parameter", "population")

# Function for Leave-one-out analysis
abc.loo <- function(empir, parami, simi, stat2rm) {
  # empiri=empir; parami=models %>% dplyr::select(-c(bottle)); simi=sim.stat; stat2rm="Dim.1_corr.time.est"
  
  abc.tmp <- abc(target = empir %>% dplyr::select(-c(stat2rm)), 
                  param = parami,
                  sumstat = simi %>% dplyr::select(-c(stat2rm)),
                  tol=0.05, 
                  transf=c("none","log"),
                  method="loclinear")
  
  # Summarize
  abc.fin <- data.table(summary(abc.tmp),
                        statrm=stat2rm)

  return(abc.fin)
}

# Go through leave-one-out
abc.out <- foreach(i=1:length(sstats), .combine = "rbind") %do% {
  # i=2
  
  print(i)
  print(sstats[i])
  
  # ABC using rejection methods
  l <- abc.loo(empir = empir, 
                 parami = models %>% dplyr::select(-c(bottle)),
                 simi = sim.stat, 
                 stat2rm = sstats[i])
}

colnames(abc.out)[1:3] <- c("data", "parameter", "population")

# Merge with best estimated with every stat
full <- data.table(data.table(tmp.out[threshold==0.05] %>% 
  dplyr::select(c(data, parameter, population)),
  statrm="standard") %>% 
  full_join(abc.out) %>% 
  pivot_wider(names_from = data, values_from = population) %>% 
  rename('Perc0.25'='Weighted 2.5 % Perc.:',
         'Perc0.975'='Weighted 97.5 % Perc.:',
         'Mean'='Weighted Mean:',
         "Mode"="Weighted Mode:",
         "Max"="Max.:",
         "Median"="Weighted Median:",
         "Min"="Min.:"))
  
### Plotting data ###

# Plot within vs yearly FST
fst.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("1.within_FST",
                        "2.Overwinter_FST",
                        "3.Multi-Year_FST")) %>%
    mutate(Year=case_when(name=="1.within_FST" ~ 1,
                          name=="2.Overwinter_FST" ~ 2,
                          name=="3.Multi-Year_FST" ~ 3)) %>% 
    mutate(bot.class=case_when(bot>=0.90 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    group_by(model, nMax, nMin, Year, bot.class) %>% 
    summarize(median.fst=median(value)) %>% 
    #filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x=as.factor(Year),
               y=median.fst, 
               fill=bot.class)) + 
    geom_boxplot() +
    geom_point(aes(x=1, y=0.002060831), color="blue", size=4) +
    geom_point(aes(x=2, y=0.00328538), color="blue", size=4) +
    geom_point(aes(x=3, y=0.003746722), color="blue", size=4) +
    theme_bw() +
    scale_y_log10() +
    labs(x="Number of Winters",
         y="log10(Median Fst)",
         fill="Bottleneck %") +
    scale_fill_brewer(palette="BuPu") +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

# Plot PCA 
pca.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("Dim.1_corr.time.est",
                        "LD1_corr.time.est",
                        "LD2_corr.time.est")) %>%
    mutate(fact=case_when(name=="LD1_corr.time.est" ~ "LD 1",
                          name=="LD2_corr.time.est" ~ "LD 2",
                          name=="Dim.1_corr.time.est" ~ "PC 1",)) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    group_by(model, nMax, nMin, fact, bot.class) %>% 
    summarize(median.pca=median(value)) %>% 
    #filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x=fact,
               y=median.pca, 
               fill=bot.class)) + 
    geom_boxplot() +
    geom_point(aes(x="LD 1", y=0.7774567), color="blue", size=4) +
    geom_point(aes(x="LD 2", y=0.2350853), color="blue", size=4) +
    geom_point(aes(x="PC 1", y=0.3686491), color="blue", size=4) +
    theme_bw() +
    labs(x="Statistic",
         y="Median r^2 Factor ~ Year",
         fill="Bottleneck %") +
    scale_fill_brewer(palette="BuPu") +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

### Leave-one-out analysis ###

# Create a data frame for the highlighted area
highlight_data <- data.table(ymin = full[statrm=="standard"]$Perc0.25,
                             ymax = full[statrm=="standard"]$Perc0.975,
                             Median = full[statrm=="standard"]$Median,
                             statrm=0,
                             parameter=c(full[statrm=="standard"]$parameter))

# LOO fig
loo <- {
  full[!statrm=="standard"] %>% 
  ggplot(.,aes(x=statrm,
               y=Median,
               ymin=Perc0.25,
               ymax=Perc0.975,
               color=parameter)) +
  geom_rect(data = highlight_data, 
          aes(xmin = -Inf, xmax = Inf,
              ymin=ymin, ymax=ymax), fill = "yellow", alpha = 0.2) +
  geom_hline(data = highlight_data, aes(yintercept=Median), linetype=2) +
  geom_linerange(size=1) +
  geom_point(size=3) +
  facet_wrap(~parameter, nrow = 2, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("#9966FF","#FF9933")) +
  labs(x="Statistic dropped",
       y="Estimated population size",
       color="") +
  theme(strip.text = element_text(face="bold", size=18),
        title = element_text(face="bold.italic", size=20),
        legend.text = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.title = element_text(face="bold", size=20))
}

# Output final plots
ggsave(fst.plot, width = 12, height = 8, file = "fst.plot.empir.new.pdf")
ggsave(pca.plot, width = 12, height = 8, file = "pca.plot.empir.new.pdf")
ggsave(loo, width = 12, height = 8, file = "loo.plot.empir.new.pdf")

### In-text information regarding bottleneck size and harmonic mean ###

# Show ABC Results for 5% threshold
tmp.out %>% 
  as.data.frame() %>%
  filter(threshold %in% 0.05)

# 5% threshold calculation
nMax=27584.1089
nMin=282.7212

# Harmonic mean
a <- c(rep(nMax, 15), rep(nMin, 2))
1/mean(1/a)

# Bottleneck
((nMax-nMin)/nMax)*100
