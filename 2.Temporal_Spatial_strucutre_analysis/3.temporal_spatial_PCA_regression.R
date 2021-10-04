### R Script for talk and paper
### 
## This script makes figure 3 of the paper

rm(list = ls())
# Load packages

#args = commandArgs(trailingOnly=TRUE)

## the 2 arguments
## Argument 1 is the Rdata with the genotype matrix
## Argument 2 is the city to analyze
## Argument 3 is the name of the analysis

#data_in=args[1]
data_in="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.ECfiltered.Rdata"

#city_target=unlist(strsplit(args[2], ",") )

#name=args[3]
name="global_set"

#load packages
library(tidyverse)
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
library(scales)
library(gmodels)

##Import the inversion mapping SNPs
inversions = read.table(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/inversion_makers_kapun.txt",
  head = T)

names(inversions)[2:3] = c("chr","pos")

## Load the genotype matrix
load(data_in)

#### Begin Measure inversions
### Measure inversion proportions
inversions %>% head

dat_filtered_t %>% 
  separate(SNP_id, 
           into = c("chr","pos"), 
           sep = "_", 
           remove = F) ->
  dat_filtered_t_chr_pos

dat_filtered_t_chr_pos$pos = as.numeric(dat_filtered_t_chr_pos$pos)

left_join(inversions, 
          dat_filtered_t_chr_pos, 
          by = c("chr","pos")) %>%
  .[complete.cases(.),] %>%
  melt(id = c("inversion", 
              "chr", 
              "pos", 
              "SNP_id", 
              "allele")) %>% 
  group_by(inversion,variable) %>%
  summarize(Inv_freq = mean(value)) %>% 
  dcast(variable~inversion, value.var = "Inv_freq") ->
  Inversion_frequencies

names(Inversion_frequencies)[1] = "sampleId" 

left_join(filtered_samps_for_analysis, Inversion_frequencies) ->
  filtered_samps_for_analysis
########## <--- End measure inverions

###
###
###

dat_filtered_t %>% ###<-- this object was created above
  .[,-which(names(.) %in% c("SNP_id"))] %>%
  t() ->
  dat_AF 

filtered_samps_for_analysis$city = gsub("Charlotttesville","Charlottesville", 
                                        filtered_samps_for_analysis$city )

filtered_samps_for_analysis$city = gsub("Odesa","Odessa", 
                                        filtered_samps_for_analysis$city )

##filtered_samps_for_analysis %>%
##  .[which(.$city %in% city_target ),] %>%
##  .$sampleId -> samps_target
##
##if("Charlottesville" %in% city_target) {
##  samps_target[-which(samps_target
##                      %in% c("VA_ch_12_spring", 
##                             "VA_ch_12_fall" ))] ->
##    samps_target
##}

dat_AF %>%
  t() %>% 
  as.data.frame -> dat_AF_samps_target

## Some characterizations of AFs and subsequent filtering
MeanAF=c()
MinAF=c()

apply(dat_AF_samps_target,
      1, FUN=mean, na.rm=TRUE ) -> MeanAF
data.frame(SNP_id = dat_filtered_t$SNP_id, MeanAF) -> MeanAF

apply(dat_AF_samps_target,
      1, FUN=min, na.rm=TRUE ) -> MinAF
data.frame(SNP_id = dat_filtered_t$SNP_id, MinAF) -> MinAF

cbind(dat_AF_samps_target, MeanAF, MinAF[-1]) -> dat_AF_samps_target

##
dat_AF_samps_target %>%
  .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MinAF > 0.001),] ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  summarise(N=n())

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  group_by(chr) %>%
  summarise(N=n())

# Remove X
dat_AF_samps_target_filtered %>%
  .[grep("X", .$SNP_id, invert = T ),] %>% 
  .[,-which(names(.) %in% c("SNP_id", "MeanAF", "MinAF"))] %>%
  t() ->
  dat_for_Analysis

colnames(dat_for_Analysis) = dat_AF_samps_target_filtered$SNP_id[grep("X", dat_AF_samps_target_filtered$SNP_id, invert = T )]


### Run up to here to begin
#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 
#### here i am making a pca object with all populations

dat_for_Analysis %>%
  as.data.frame() %>% 
  PCA(scale.unit = F, graph = F, ncp = 100) ->
  PCA_object

save(PCA_object, file = "PCA.object.all.Rdata")
load("./PCA.object.all.Rdata")

#ectracts variace
PCA_object$eig[1,2] -> D1VE
PCA_object$eig[2,2] -> D2VE

#Plot PCA
## Plot object
print("115")
PCA_object$ind$coord[,c(1,2)] %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.),
         Run = name) %>% 
  left_join(., filtered_samps_for_analysis ) -> PCA_table

save(PCA_table, file = "PCA_table.allchroms.Rdata")
load("./PCA_table.allchroms.Rdata")

####
#### Make the centroid of each population
label.df_PCA <- PCA_table %>% 
  group_by(city) %>% 
  summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2) ) 

colorset = c(brewer.pal(n = 12, name = "Set3"),
             "darkgoldenrod1",
             "brown2",
             "hotpink4",
             "olivedrab4"
)

### Draw the plot
ggplot() + 
  geom_point(data  =PCA_table,
             shape=21,
             aes(x=Dim.1, y=Dim.2, 
                 fill = city,
                 #shape = continent,
                 size = 3,  
                 alpha = 0.8))  + 
  geom_text_repel(data = label.df_PCA, size = 3 ,  
                  box.padding = 0.5, max.overlaps = Inf,
                  alpha = 1, 
                  aes(x=Dim.1, 
                      y=Dim.2, 
                      label = city))  + 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values =  colorset) +
  xlab(paste("PC1 (", round(D1VE, 2), "% VE)" , sep = "")) +
  ylab(paste("PC2 (", round(D2VE, 2), "% VE)", sep = "")) ->
  #scale_shape_manual(values = 
  #                     c(21,22,23,24)) -> 
  PCA_12 #<<<<<< Panel A


#Save the plot
ggsave(PCA_12, 
       file = paste(name, "PCA", "pdf" , sep = ".") ,
       width = 4,
       height = 4)

############### Panel B
############### 

PCA_table %>%
  #.[,c("Dim.1","Dim.2", 
  #     "country","city","collectionDate","year"
  #    )] %>%
  #melt(id = c("country","city","collectionDate","year")) %>%
  ggplot(aes(
    x=  as.Date(collectionDate, 
                format = c("%m/%d/%Y")),
    y=Dim.1,
    fill = city,
  )) +
  geom_smooth(method = "lm", 
              se = F,
              size = 1,
              aes(color = city)) +
  geom_point(
    size = 2,  
    alpha = 0.8,
    shape = 21)  + 
  theme_bw() +
  facet_wrap(~city, 
             scales = "free", 
             nrow =3
             ) +
  theme(legend.position = "none")  +
  scale_fill_manual(values =  colorset) +
  scale_color_manual(values =  colorset) +
  ylab("Global PC 1") +
  xlab("Collection Year") + 
  scale_x_date(labels = date_format("%Y"), date_breaks = "1 year") +
  scale_shape_manual(values = c(21,22)) -> 
  PC1_v_year

ggsave(PC1_v_year, 
       file = paste(name, "pc1_v_date", "pdf" , sep = ".") ,
       width = 8,
       height = 4)

#######
#######

load("./PCA_table.allchroms.Rdata")

PCA_table$city %>% unique() -> city_guide

cor_list = list()

for(i in 1:length(city_guide)){
  
  tmp = cor.test( ~ Dim.1 + year, data = PCA_table[which(PCA_table$city == city_guide[i]),] )
  tmp_df = data.frame(pop = city_guide[i],
             cor = tmp$estimate,
             p.val= tmp$p.value)
  
  cor_list[[i]] = tmp_df
  
}

cor_df = do.call(rbind, cor_list)

cor_df[order(cor_df$pop),]

cor_df %>%
  .[which(.$p.val < 0.1),] %>%
  .[order(.$pop),] %>%
  mutate(R2 = cor^2)


cor_df %>%
  .[which(.$p.val < 0.1),] %>%
  .$pop %>%
  unique() -> cities_sigP
  

#PCA eucledian distance analysis
######################
###################### Shuffle analysis
cities_sigP = c("Munich", "Linvilla", "Yesiloz", "Odessa",
                "Yalta", "Charlottesville" )

PCA_table %>% 
  .[which(.$city %in% cities_sigP),] %>%
  group_by(city, year) %>% 
  summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2)) -> mean_pca_dist

shuffle_list = list()

for(i in 1:length(cities_sigP)){

  print(cities_sigP[i])
  mean_pca_dist %>%
    .[which(.$city == cities_sigP[i]),] ->
    data_in
  
    data_in %>%
    .[,-c(1,2)] %>%
    as.matrix() %>%
    spDists( longlat=F) %>%
    as.matrix() -> mean_pca_matrix

  
    data_in %<>%
    .[which(.$city == cities_sigP[i]),] %>%
    mutate(name_year = paste(year, sep = "_"))
  
  dimnames(mean_pca_matrix) <- list(data_in$name_year, 
                                    data_in$name_year) 
  
  xy <- t(combn(colnames(mean_pca_matrix), 2))
  
  data.frame(xy, dist=mean_pca_matrix[xy]) %>%
    mutate(type = paste(name, "actual", sep = "."),
           city = cities_sigP[i])  -> 
    distance_output.actual
  
  random_out =list()
  
  for(k in 1:1000){
    
    print(paste(cities_sigP[i], k, sep = " "))
    
    data_in %>% 
      mutate(
        random_y=data_in$year[shuffle(data_in$year)]  
      ) %>%
      group_by(random_y) %>% 
      summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2)) -> random_shuffle
    
    random_shuffle %>%
      .[,-c(1,2)] %>%
      as.matrix() %>%
      spDists( longlat=F) %>%
      as.matrix() -> mean_random_shuffle
    
    random_shuffle %<>%
      mutate(name_year = paste(random_y, sep = "_"))
    
    dimnames(mean_random_shuffle) <- list(random_shuffle$name_year, 
                                          random_shuffle$name_year) 
    
    xyr <- t(combn(colnames(mean_random_shuffle), 2))
    
    data.frame(xyr, dist=mean_random_shuffle[xyr]) %>%
      mutate(type = paste(name, "random", i, sep = "."),
             city = cities_sigP[i]) ->
      distance_output.random
    
    random_out[[k]] = distance_output.random
    
  }
  
  random_df.shuffle = do.call(rbind, random_out)
  
  shuffle_list[[i]] = rbind(distance_output.actual, random_df.shuffle)
  
}

shuffle_df = do.call(rbind, shuffle_list)

### estimate means
### 

shuffle_df$type[grep("random", shuffle_df$type)] = "all.random"

## run MannU.test
t_list = list()

for(i in 1:length(cities_sigP)){
  
  tmp = wilcox.test( dist ~ type, 
                data = shuffle_df[which(shuffle_df$city == cities_sigP[i]),] )
  
  tmp_df = data.frame(pop = cities_sigP[i],
                      p.val= tmp$p.value)
  
  t_list[[i]] = tmp_df
  
}

t_df = do.call(rbind, t_list)
