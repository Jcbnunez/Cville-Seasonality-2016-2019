# Compiles and analyze genomic data from SLiM and PLINK
# Connor Murray 10.13.2021
# ijob -A berglandlab_standard --mem=10G -p standard -c 4
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(forcats)
library(viridis)
library(gganimate)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)

# VCF output folder
setwd("/scratch/csm6hg/slim_bottleneck/")

# Load data by chromosome
load("Cville_2L.ECfiltered.Rdata")
L2 <- o
load("Cville_2R.ECfiltered.Rdata")
R2 <- o 
load("Cville_3L.ECfiltered.Rdata")
L3 <- o
load("Cville_3R.ECfiltered.Rdata")
R3 <- o

# Merge all chromosomes
test <- merge(merge(L2, R2, by="row.names"), 
              merge(L3, R3, by="row.names"), 
              by="row.names")

# QC
dim(test)

# Convert to data table
af <- data.table(test)
samps <- unique(data.table(af %>% filter(Row.names.x %in% filtered_samps_for_analysis[MeanEC >30]$sampleId))$Row.names.x)

# Transpose
af %>%
  filter(Row.names.x %in% filtered_samps_for_analysis[MeanEC >30]$sampleId) %>% 
  select(-c(colnames(af)[which(colnames(af) %like% "names")])) %>%
  t() %>% 
  as.data.frame -> dat_AF_samps_target

# Filtering by AF and missingness 
MeanAF=c()
MinAF=c()

# Mean AF
apply(dat_AF_samps_target, 1, FUN=mean, na.rm=TRUE) -> MeanAF
data.frame(SNP_id = 1:length(MeanAF), MeanAF) -> MeanAF

# Minimum AF
apply(dat_AF_samps_target, 1, FUN=min, na.rm=TRUE) -> MinAF
data.frame(SNP_id = 1:length(MinAF), MinAF) -> MinAF

cbind(dat_AF_samps_target, MeanAF, MinAF[-1]) -> dat_AF_samps_target

# Filtering Polymorphic
dat_AF_samps_target %>%
  .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MinAF > 0.001),] ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

dim(dat_AF_samps_target)
dim(dat_AF_samps_target_filtered)

# Missingness
count_NA = function(x){
  return(sum(is.na(x)))
}
MissDat=c()

pool_cols = grep("SNP_id|MeanAF|MinAF", colnames(dat_AF_samps_target_filtered), invert = T)

apply(dat_AF_samps_target_filtered[,pool_cols], 1, FUN=count_NA ) -> MissDat

n_pools = length(pool_cols)

data.frame(SNP_id = dat_AF_samps_target_filtered$SNP_id, missing_rate = c(MissDat/n_pools) ) -> MissDat

cbind(dat_AF_samps_target_filtered, MissDat[-1]) -> dat_AF_samps_target_filtered

# This samples only polymorphic sites
dat_AF_samps_target_filtered %>%
  filter(missing_rate < 0.01) -> af 

# Add sample names
colnames(af)[1:31] <- samps

# Run PCA
af %>%
  select(-c(colnames(af)[32:35])) %>%
  t() %>% 
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  pca

# Include sample metadata
pca.dt <- data.table(pca$ind$coord, 
                     Sample=colnames(t(pca$ind$dist))) %>% 
          left_join(filtered_samps_for_analysis %>% 
                    select(Sample=sampleId, Year=year, population=city, yday, MeanEC)) %>% 
  mutate(dd.sim = case_when(yday < 200 ~ "1.within",
                            yday >= 200 ~ "2.overwinter"))

# Add correlation coefficients
pca.dt <- data.table(pca.dt %>% 
           mutate(r2pc1 = summary(lm(data = pca.dt, Dim.1 ~ Year))$r.squared,
                  r2pc2 = summary(lm(data = pca.dt, Dim.2 ~ Year))$r.squared,
                  r2pc3 = summary(lm(data = pca.dt, Dim.3 ~ Year))$r.squared))

pdf("pca.charlottesville.include2L.filtered.pdf")

# Plot PCA
ggplot(pca.dt) + 
  geom_point(aes(x = Dim.1, 
                 y = Dim.2, 
                 color = as.factor(Year)),
             size=8) +
  labs(x=paste("PC1 ", "(", round((pca$eig[1,2]), digits=2), "%)", sep=""), 
       y=paste("PC2 ", "(", round((pca$eig[2,2]), digits=2), "%)", sep=""),
       title="Including 2L",
       color="Year") + 
  theme_bw() +
  theme(title = element_text(face="bold", size=20),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()

# Save output
save(pca.dt, file = "PCA.object.include2L.filtered.Rdata")

### EXCLUDE 2L ###

# Make data.table excluding 2L
af2 <- data.table(af[rownames(af) %like% "2L_" == FALSE,])

# QC
dim(af2)

# Run PCA without 2L
af2 %>% 
  select(-c(colnames(af2)[32:35])) %>%
  t() %>% 
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  pca2

# Include sample metadata without 2L
pca.dt2 <- data.table(pca2$ind$coord, 
                      Sample=colnames(t(pca$ind$dist))) %>% 
  left_join(filtered_samps_for_analysis %>% 
              select(Sample=sampleId, Year=year, population=city, yday, MeanEC)) %>% 
  mutate(dd.sim = case_when(yday < 200 ~ "1.within",
                            yday >= 200 ~ "2.overwinter"))

# Add correlation coefficients
pca.dt2 <- data.table(pca.dt2 %>% 
            mutate(r2pc1 = summary(lm(data = pca.dt2, Dim.1 ~ Year))$r.squared,
                   r2pc2 = summary(lm(data = pca.dt2, Dim.2 ~ Year))$r.squared,
                   r2pc3 = summary(lm(data = pca.dt2, Dim.3 ~ Year))$r.squared))

pdf("pca.charlottesville.exclude2L.filtered.pdf")

# Plot PCA
ggplot(pca.dt2) + 
  geom_point(aes(x = Dim.1, 
                 y = Dim.2, 
                 color = as.factor(Year)),
             size=8) +
  labs(x=paste("PC1 ", "(", round((pca2$eig[1,2]), digits=2), "%)", sep=""), 
       y=paste("PC2 ", "(", round((pca2$eig[2,2]), digits=2), "%)", sep=""),
       title="Excluding 2L",
       color="Year") + 
  theme_bw() +
  theme(title = element_text(face="bold", size=20),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=18),
        legend.background = element_blank(),
        strip.text =element_text(face="bold", size=18),
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20))

dev.off()

# Save output
save(pca.dt2, file = "PCA.object.exclude2L.filtered.Rdata")