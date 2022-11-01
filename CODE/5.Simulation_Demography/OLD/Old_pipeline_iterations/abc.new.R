# Libraries
library(abc)
library(abc.data)
library(data.table)
library(foreach)
library(tidyverse)
library(ggrepel)

# Working directory
setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/slim-abc/Data/")

# Empirical diversity data
data <- data.table(na.omit(read.csv("MapJune2020.recode.finalsnpset.F1Wilds.diversity.csv", 
                                    header = T)))

data <- data.table(data %>% select(tW.stan.pop.1, tP.stan.pop.1, tajima.D.pop.1, SC, chr) %>% 
                     mutate(tP.stan=as.numeric(tP.stan.pop.1),
                            tW.stan=as.numeric(as.character(tW.stan.pop.1)),
                            tajimaD=as.numeric(tajima.D.pop.1),
                            SC=SC))

# Plot data
ggplot(data[!SC=="OO"], aes(x=SC, y=tP.stan, color=as.factor(SC))) + geom_boxplot()

# Mean of clonal lineages
# data <- data.table(data %>% group_by(SC) %>% select(tP.stan, tW.stan, tajimaD) %>% summarise_each(funs = mean))

# Mean and CI of clonal lineages
data <- data[,list(tP.stan=mean(tP.stan, na.rm=T),
                   tW.stan=mean(tW.stan, na.rm=T),
                   tajimaD=mean(tajimaD, na.rm=T),
                   tP.stan.lci=quantile(tP.stan, 0.025, na.rm=T), 
                   tP.stan.uci=quantile(tP.stan, 0.975, na.rm=T),
                   tW.stan.lci=quantile(tW.stan, 0.025, na.rm=T), 
                   tW.stan.uci=quantile(tW.stan, 0.975, na.rm=T),
                   tajimaD.lci=quantile(tajimaD, 0.025, na.rm=T), 
                   tajimaD.uci=quantile(tajimaD, 0.975, na.rm=T)),
                          list(SC)]

# Old diversity
data2 <- data.table(readRDS("filtered-10rep-15samp.rds"))
colnames(data2)[1:6] <- c("tW", "tW.stan", "tP", "tP.stan", "tP.eff", "tajimaD")

# Plot data
ggplot(data2[!SC=="OO"], aes(x=SC, y=tP.stan, color=as.factor(SC))) + geom_point()

data2 <- data.table(data2 %>% group_by(SC) %>% select(tP.stan, tW.stan, tajimaD)%>% 
                     summarise_each(funs = mean))

# Simulation diversity data
sim <- data.table(rbind(read.csv("constant.output.10.20.csv")))

colnames(sim)[1:6] <- c("tW", "tW.stan", "tP", "tP.stan", "tP.eff", "tajimaD")
sim[, seed:=tstrsplit(sim$vcf, "_")[10]]

# Plot data
ggplot(sim[EG < 1.4], aes(x=gen, y=tP.stan, color=as.factor(K))) + geom_boxplot()

# Model information
sim[, model:=paste("constant", EG, K, gen, sep="_")]
models <- as.character(sim$model)

# Sim stats
sim.stat <- data.table(sim %>% select(tP.stan, tW.stan, tajimaD))

# Foreach abc stats
out <- foreach(i=1:length(data$SC), .combine = "rbind") %do% {
  
  # ABC using rejection methods
  tmp <- abc(target = data[i,-c(1, 5:10)],
            param=sim[,c(24)], 
            sumstat=sim.stat,
            tol=0.1, 
            method="rejection")
  
  # Compile data
  tmp <- data.table(summary(tmp), 
                    SC=data[i]$SC,
                    data[i, c(2:10)])
  
}

# Fix data
colnames(out)[1:4] <- c("data", "sim.gen", "Gen", "SC")
out[data %in% c("2.5% Perc.:","97.5% Perc.:","Mean:")]

out <- data.table(pivot_wider(out, names_from = data, values_from = Gen))
colnames(out)[12:18] <- c("Min.", "Perc0.25", "Median", "Mean", "Mode", "Perc0.975", "Max.") 

# Age and diversity plot
ggplot(out, aes(x=Mean, y=log10(tP.stan), color=SC, label=SC)) +
  geom_pointrange(aes(ymin=log10(tP.stan.uci), ymax=log10(tP.stan.lci)), fatten = .5, size = 1.2) +
  geom_pointrange(aes(xmin=Perc0.25, xmax=Perc0.975), fatten = .5, size = 1.2) +
  geom_label(size=10, alpha=0.8) +
  labs(x="Age (Generations)", y="log10(Nucleotide diversity)") +
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") 
