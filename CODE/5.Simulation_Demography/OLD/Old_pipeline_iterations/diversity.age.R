# Nucleotide diversity from VCF
# 10.20.2020

# Libraries
library(PopGenome)
library(data.table)
library(foreach)
library(tidyverse)
library(viridis)
library(cowplot)
library(ggforce)

# Empirical diversity data
data <- data.table(read.csv("C:/Users/Conno/Desktop/Fall2020/Berglandlab/slim-abc/Data/div.fig.csv", header = T))
data.w <- data %>% pivot_wider(names_from = factor, values_from = data)

# Diversity and simulation
ggplot(data=data.w, aes(x=mean, y=mean.pi, color=SC, label=SC)) +
   geom_pointrange(aes(ymin=lci.pi, ymax=uci.pi),
                   fatten = .5, size = 1.2) +
   geom_pointrange(aes(xmin=lci, xmax=uci),
                   fatten = .5, size = 1.2) +
   geom_label(size=14) +
   labs(x="Age (Generations)", y="Nucleotide diversity",
        color="Clonal lineage") +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=16),
         axis.text.y = element_text(face="bold", size=16),
         axis.title.x = element_text(face="bold", size=18),
         axis.title.y = element_text(face="bold", size=18),
         legend.background = element_blank(),
         legend.text = element_text(size=14),
         legend.title = element_text(face="bold", size=16),
         legend.position = "none") 

# Simulation diversity data
sim <- data.table(read.csv("C:/Users/Conno/Desktop/Fall2020/Berglandlab/slim-abc/Data/constant.test.output.csv"))

colnames(sim)[1:6] <- c("tW", "tW.stan", "tP", "tP.stan", "tP.eff", "tajimaD")
sim[, seed:=tstrsplit(sim$vcf, "_")[10]]

# Diversity long format
nuc.div <- na.omit(sim[,list(avg.pi = mean(tP.stan),
                                avg.watterson = mean(tW.stan),
                                avg.tajima = mean(tajimaD)),
                           list(EG, K, gen)])

sort(unique(sim$gen))
sort(unique(sim$K))

# Labels for heatmaps
label <- data.table(x=c(38, 9, 47, 79),
                    y=c(14, 15, 14, 14),
                    label=c("A", "B", "C", "D"))

# Theta pi
a <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.pi)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=2.1e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=5.3e-7), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=2.5e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=4.4e-6), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Exponential growth", y="K", fill="Theta pi") +
   scale_fill_viridis(option = "A") +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_blank(),
         axis.title.y = element_text(face="bold", size=10),
         axis.ticks.x = element_blank()) 

# Theta Watterson
b <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.watterson)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=7.1e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=1.7e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=8.3e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=1.3e-5), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Exponential growth", y="K", fill="Theta Wa") +
   scale_fill_viridis(option = "A") +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_blank(),
         axis.title.y = element_text(face="bold", size=10),
         axis.ticks.x = element_blank())

# Tajima's D
c <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.tajima)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=-2.70), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=-2.55), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=-2.66), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=-2.59), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Generation", y="K", fill="Tajima's D") +
   scale_fill_viridis(option = "A") +
   scale_x_discrete(breaks=c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)) +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=8, vjust = 0.4),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_text(face="bold", size=10),
         axis.title.y = element_text(face="bold", size=10)) 

plot_grid(a,b,c, nrow=3)

