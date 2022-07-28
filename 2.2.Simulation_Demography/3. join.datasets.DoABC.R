### 3. Join analyses for ABC
### 

library(foreach)
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(viridis)
library(gganimate)
library(abc)
library(abc.data)
library(lme4)
library(ggrepel)
library(ggforce)
library(cowplot)

####
####
root <- "/scratch/yey2sn/Overwintering_ms/20.Connor_sims/Parsed_sims"
files <- system( paste("ls ", root), intern = T )

joint.dat <- foreach(i=1:length(files), .combine = "rbind", .errorhandling = "remove" )%do%{
 
  message(paste(i, "/",length(files)))
  tmp <- get(load( paste( "Parsed_sims/" , files[i], sep = "" ) ))
  
}

#####

# Load empirical data: obs_values.final
load("/project/berglandlab/connor/bottleneck/obs_values.final.Rdata")

# Change names of empirical data
empir <- data.table(obs_values.final %>% 
                      pivot_wider(names_from="Obs_var", values_from="Median_Val"))
empir <- empir[,c(1:8,11,19,25,26,32:34)]
colnames(empir) <- as.character(tstrsplit(colnames(empir),"=")[[2]])


### work with sim data
mean.sim <- joint.dat
mean.sim <- mean.sim[,c(1:4, 15:17,30,32,12:14,6, 61:63,49:51)][nMax>1000]

colnames(mean.sim)[8:19] <- c("LD1_corr.time.est", "LD2_corr.time.est", 
                              "comp_1_percentage_of_variance", "comp_2_percentage_of_variance",
                              "comp_3_percentage_of_variance","50%_AF_var",
                              "0_Median.standarized.PCAEucDist","1_Median.standarized.PCAEucDist", 
                              "2_Median.standarized.PCAEucDist",
                              "3.Multi-Year_FST", "2.Overwinter_FST", "1.within_FST")


# Model information
models <- data.table(nMax=as.numeric(tstrsplit(mean.sim$model, "_")[[1]]),
                     nMin=as.numeric(tstrsplit(mean.sim$model, "_")[[2]])) %>% 
  mutate(bottle=((nMax-nMin)/nMax)*100)

# Statistics to use
sstats <-  c(colnames(empir))
sstats <- sstats[-c(6:8,10,11,13)]

# Empirical stats
empir <- empir %>% dplyr::select(sstats)

# ABC stats
sim.stat <- mean.sim %>% dplyr::select(sstats) %>% 
  relocate(colnames(empir))

# Thresholds to test in ABC
thresh <- seq(from=0.01, to = 0.15, by = 0.001)

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


######
######
###### ---> plots
######
######


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

out %>% 
  as.data.frame() %>%
  filter(threshold %in% seq(from=0.09, to = 0.11, by = 0.01))


#pdf("thresholds.abc.rejection.exclude2L.nopcdims.bottleneck.xlim.more.pdf")

# ABC and tolerance thresholds
ggplot(out[!parameter=="bottle"], 
       aes(x=threshold, 
           y=as.numeric(Mean),
           ymin=as.numeric(Perc0.25),
           ymax=as.numeric(Perc0.975),
           color=parameter)) +
  geom_pointrange(position = "dodge") +
  labs(x="Thresholds", 
       y="Bottleneck size %", 
       color="Parameter") +
  theme_classic() +
  scale_fill_viridis(option = "magma", begin=0, end=1) +
  theme(title = element_text(face="bold", size=15),
        legend.text = element_text(face="bold", size=14),
        legend.title = element_text(face="bold", size=15),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=15),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16))

#dev.off()

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
    mutate(bot.class=case_when(bot>=0.95 ~ "Above 95%",
                               bot>=0.9 & bot <0.95 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    ggplot(., 
           aes(x=as.factor(Year), 
               y=value, 
               group=interaction(replicate,model), 
               color=bot)) + 
    geom_point(size=2.2, alpha=0.6) +
    geom_line(alpha=0.1, linetype=1) +
    facet_wrap(~bot.class, nrow=1) +
    theme_bw() +
    labs(x="Number of Winters (sim. year)",
         y="Fst",
         color="Bottleneck %") +
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
    ggplot(., 
           aes(x=as.factor(fact), y=value, group=interaction(model,replicate), color=bot)) + 
    geom_point(size=1, alpha=0.4) +
    geom_line(alpha=0.1) +
    facet_wrap(~bot.class, nrow=1) +
    theme_bw() +
    labs(x="Principal component or Linear Disriminant Axis",
         y="r^2 Factor ~ Year",
         color="Bottleneck %") +
    scale_color_viridis(begin = 1, end=0) +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18, angle = -40),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

# EUC plot
euc.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("0_Median.standarized.PCAEucDist", "1_Median.standarized.PCAEucDist", 
                        "2_Median.standarized.PCAEucDist")) %>%
    mutate(fact=case_when(name=="0_Median.standarized.PCAEucDist" ~ "EUC 0",
                          name=="1_Median.standarized.PCAEucDist" ~ "EUC 1",
                          name=="2_Median.standarized.PCAEucDist" ~ "EUC 2")) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    ggplot(., 
           aes(x=as.factor(fact), y=value, group=interaction(model,replicate), color=bot)) + 
    geom_point(size=1, alpha=0.4) +
    geom_line(alpha=0.1) +
    facet_wrap(~bot.class, nrow=1) +
    theme_bw() +
    labs(x="Number of winters (Years)",
         y="PCA Euclidean Distance",
         color="Bottleneck %") +
    scale_color_viridis(begin = 1, end=0) +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18, angle = -40),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))}

# Percent explained plot
perc.plot <- {mean.sim %>% 
    mutate(bot=(nMax-nMin)/nMax) %>% 
    pivot_longer(cols=c("comp_1_percentage_of_variance", "comp_2_percentage_of_variance", 
                        "comp_3_percentage_of_variance")) %>%
    mutate(fact=case_when(name=="comp_1_percentage_of_variance" ~ "PC % 1",
                          name=="comp_2_percentage_of_variance" ~ "PC % 2",
                          name=="comp_3_percentage_of_variance" ~ "PC % 3")) %>% 
    mutate(bot.class=case_when(bot>=0.9 ~ "Above 90%",
                               bot<0.9 & bot > 0.5 ~ "50-90%",
                               bot<=0.5 & bot>0.1~ "10-50%",
                               bot<=0.1 ~ "0-10%")) %>% 
    ggplot(., 
           aes(x=as.factor(fact), y=value, group=interaction(model,replicate), color=bot)) + 
    geom_point(size=1, alpha=0.4) +
    geom_line(alpha=0.1) +
    facet_wrap(~bot.class, nrow=1) +
    theme_bw() +
    labs(x="PC Number",
         y="Percentage PC variance explained",
         color="Bottleneck %") +
    scale_color_viridis(begin = 1, end=0) +
    theme(strip.text = element_text(face="bold", size=18),
          title = element_text(face="bold.italic", size=20),
          legend.text = element_text(face="bold", size=16),
          legend.title = element_text(face="bold", size=16),
          axis.text.x = element_text(face="bold", size=18, angle = -40),
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
    ggplot(., 
           aes(x=value, group=bot.class, color=bot.class)) + 
    geom_density() +
    theme_bw() +
    labs(x="Allele frequency variance",
         y="Density",
         color="Bottleneck %") +
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

ggsave(sse.plot, file = "sse.plot.pdf")

#pdf("../../figures/composite.abc.new2.pdf", width = 35, height = 16)
plot_grid(fst.plot, pca.plot, af.plot, nrow = 1, labels = c("A.", "B.", "C."))
dev.off()

a <- c(rep(25000, 10), rep(1500,2))
1/mean(1/a)

# All frequency data
filenames1 <- list.files(path = "/project/berglandlab/connor/bottleneck/data3", pattern = ".txt$")
