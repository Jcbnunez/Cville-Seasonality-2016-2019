# ABC analysis and plotting of Drosophila simulations
# Connor S. Murray
# 10.31.2022
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
setwd("/project/berglandlab/connor/bottleneck")

# Simulation data
joint.dat <- readRDS("parsed.data.Rdata")

# Restrict to analyses 
restrict <- data.table(fread("model_paramList_fin") %>% 
                         mutate(model=paste(nMax, nMin, sep="_")))

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
#sstats <- sstats[-c(6:8,10,11,13)]

# Empirical stats
empir <- empir %>% dplyr::select(sstats)

# ABC stats
sim.stat <- mean.sim %>% dplyr::select(sstats) %>% 
  relocate(colnames(empir))

# Thresholds to test in ABC
thresh <- seq(from=0.01, to = 0.15, by = 0.01)

# Go through several thresholds
tmp.out <- foreach(i=1:length(thresh), .combine = "rbind") %do% {
  
  # ABC using rejection methods
  tmp <- abc(target = empir, 
             param = models %>% dplyr::select(-c(bottle)), 
             sumstat = sim.stat,
             tol=thresh[i], 
             method="rejection")
  
  # Compile data
  tmp <- data.table(summary(tmp), 
                    empir,
                    threshold=thresh[i])
  
  # Progress message
  print(paste("Threshold:", thresh[i], sep=" "))
  
  # Finish
  return(tmp)
}

### Plotting data ###

# Fix data
colnames(tmp.out)[1:3] <- c("data", "parameter", "population")
tmp.out[data %in% c("2.5% Perc.:","97.5% Perc.:","Mean:")]

out <- data.table(pivot_wider(tmp.out, names_from = data, values_from = population)) %>% 
  rename('Perc0.25'='2.5% Perc.:',
         'Perc0.975'='97.5% Perc.:',
         'Mean'='Mean:',
         "Mode"="Mode:",
         "Max"="Max.:",
         "Median"="Median:",
         "Min"="Min.:")

# Show ABC Results for 10% threshold
out %>% 
  as.data.frame() %>%
  filter(threshold %in% 
           seq(from=0.09, 
               to = 0.11, 
               by = 0.01))

# Scale simulation
sim.stat.sc <- data.table(sim.stat %>% 
             mutate_at(c(colnames(empir)),
                       ~(./mad(.)) %>% as.vector))

# Scale empir
mad.sc <- data.table(sim.stat %>% 
            summarize_at(c(colnames(empir)),
                        ~(mad(.)) %>% as.vector))
empir.sc <- data.table(empir/mad.sc)

# Sum of squares error generator 
sim.stat %>% 
  as.matrix(.) -> a
empir %>% 
  as.matrix(.) -> b

# Sum of squares
dt <- data.table(data.table(models,
                            sse=as.numeric(a %>% 
                              sweep(., 2, b) %>% 
                                .^2 %>%
                               rowSums(.))) %>% 
                   group_by(nMin, nMax) %>% 
                   summarise(sse.mean=mean(sse)))

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
    summarize(mean.fst=mean(value)) %>% 
    filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x=as.factor(Year),
               y=mean.fst, 
               fill=bot.class)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x="Number of Winters",
         y="Mean Fst",
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
                        "Dim.2_corr.time.est",
                        "Dim.3_corr.time.est",
                        "LD1_corr.time.est",
                        "LD2_corr.time.est")) %>%
    mutate(fact=case_when(name=="LD1_corr.time.est" ~ "LD 1",
                          name=="LD2_corr.time.est" ~ "LD 2",
                          name=="Dim.1_corr.time.est" ~ "PC 1",
                          name=="Dim.2_corr.time.est" ~ "PC 2",
                          name=="Dim.3_corr.time.est" ~ "PC 3")) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    group_by(model, nMax, nMin, fact, bot.class) %>% 
    summarize(mean.pca=mean(value)) %>% 
    filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x=fact,
               y=mean.pca, 
               fill=bot.class)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x="Statistic",
         y="r^2 Factor ~ Year",
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

# EUC plot
euc.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("0_Median.standarized.PCAEucDist", 
                        "1_Median.standarized.PCAEucDist", 
                        "2_Median.standarized.PCAEucDist")) %>%
    mutate(fact=case_when(name=="0_Median.standarized.PCAEucDist" ~ "EUC 0",
                          name=="1_Median.standarized.PCAEucDist" ~ "EUC 1",
                          name=="2_Median.standarized.PCAEucDist" ~ "EUC 2")) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    group_by(model, nMax, nMin, fact, bot.class) %>% 
    summarize(mean.euc=mean(value)) %>% 
    filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x=fact,
               y=mean.euc, 
               fill=bot.class)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x="Statistic",
         y="PCA Euclidean Distance",
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

# Percent explained plot
af.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("50%_AF_var")) %>%
    mutate(fact=case_when(name=="50%_AF_var" ~ "AF Var")) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    group_by(model, nMax, nMin, fact, bot.class) %>% 
    summarize(mean.afvar=mean(value)) %>%
    filter(nMin>=1000) %>% 
    ggplot(., 
           aes(x="",
               y=mean.afvar, 
               fill=bot.class)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x="",
         y="Median Allele Frequency Variance",
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

# Plot Sum of squares error
sse.plot <- {ggplot() +
    geom_raster(data=dt,
                aes(x=as.factor(nMin),
                    y=as.factor(nMax),
                    fill=sse.mean)) +
    labs(x="Minimum population size", 
         y="Maximum population size", 
         fill="Mean SSE") +
    theme_classic() +
    scale_fill_viridis(option = "magma", begin=0, end=1) +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=12),
          legend.title = element_text(face="bold", size=15),
          legend.background = element_blank(),
          strip.text =element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=18, angle = -40),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))}

dt.ma <- data.table(dt %>% 
            mutate(model=paste(nMax, nMin, sep="_"),
                   bot=((nMax-nMin)/nMax)*100,
                   bot.class=case_when(bot>=90 ~ "Above 90%",
                                       bot<90 & bot>50 ~ "50-90%",
                                       bot<=50 & bot>10 ~ "10-50%",
                                       bot<=10 ~ "0-10%")))

# Sum of Squares error matrix
sse.plot2 <- {dt.ma[bot>0] %>% 
    ggplot(.,
        aes(x=bot,
            y=sse.mean,
            color=nMin)) +
    geom_point(size=2) +
    geom_smooth(method = "gam", size=2, color="steelblue4") +
    labs(x="Bottleneck size (%)", 
         y="Mean Sum of Squares Error", 
         color="Minimum population size") +
    theme_classic() +
    scale_color_viridis(option = "magma", begin=0, end=1) +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=12),
          legend.title = element_text(face="bold", size=15),
          legend.background = element_blank(),
          strip.text =element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))}

# SSE histogram
sse.plot3 <- {dt.ma[bot>0] %>% 
    ggplot(.,
           aes(x=sse.mean,
               fill=bot.class)) +
    geom_histogram() +
    labs(x="Mean Sum of Squares Error",
         y="Number of Simulations",
         fill="Bottleneck size") +
    theme_classic() +
    scale_fill_brewer(palette = "BuPu") +
    theme(title = element_text(face="bold", size=15),
          legend.text = element_text(face="bold", size=12),
          legend.title = element_text(face="bold", size=15),
          legend.background = element_blank(),
          strip.text =element_text(face="bold", size=15),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20))}

# Save Plots
ggsave(sse.plot, width = 8, height = 8, file = "sse.plot.new.pdf")
ggsave(sse.plot2, width = 10, height = 8, file = "sse.plot.reg.new.pdf")
ggsave(sse.plot3, width = 10, height = 8, file = "sse.plot.hist.new.pdf")

pdf("composite.abc.new3.pdf", width = 35, height = 16)
plot_grid(fst.plot, pca.plot, euc.plot, nrow = 1, labels = c("A.", "B.", "C."))
dev.off()

pdf("composite.fst.pdf", width=10, height=8)
fst.plot
dev.off()

pdf("composite.pca.pdf", width=10, height=8)
pca.plot
dev.off()

pdf("composite.euc.pdf", width=10, height=8)
euc.plot
dev.off()

pdf("composite.af.pdf", width=10, height=8)
af.plot
dev.off()

# Harmonic mean calculation
a <- c(rep(16043.8, 15), rep(1084.6,2))
1/mean(1/a)

a <- c(rep(10000, 15), rep(3000,2))
1/mean(1/a)

(16043.8-1084.4)/16043.8
(10000-3000)/10000
