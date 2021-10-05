### Run Analysis of PCA
### 
rm(list = ls())

#load data
# This R object was premade in script 1 ==> 1.Import_GDStoR.r
data_in="/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/PCA.object.all.Rdata"
inmeta="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata"

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
library(data.table)

### load data
load(data_in)
load(inmeta)

samps = samps_EFFCOV

samps$city = gsub("Charlotttesville","Charlottesville", samps$city)
samps$city = gsub("Odesa","Odessa", samps$city )
samps$city[grep("Yesiloz", samps$city )] = "Yesiloz"
samps$city[grep("Chornobyl", samps$city )] = "Chernobyl"
samps$city[grep("Kyiv", samps$city )] = "Kyiv"

#ectracts variace
PCA_object$eig[1,2] -> D1VE
PCA_object$eig[2,2] -> D2VE
PCA_object$eig[3,2] -> D3VE

#Plot PCA
## Plot object
PCA_object$ind$coord[,c(1:20)] %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>% 
  left_join(., samps ) -> PCA_table

PCA_table$collectionDate = as.Date(PCA_table$collectionDate, format = "%m/%d/%Y")
####
#### Make the centroid of each population
label.df_PCA <- PCA_table %>% 
  group_by(city) %>% 
  summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2) ) 

#colorset = c(brewer.pal(n = 12, name = "Set3"),
#             "darkgoldenrod1",
#             "brown2",
#             "hotpink4",
#             "olivedrab4"
#)

### Draw the plot
ggplot() + 
  geom_point(data  =PCA_table,
             shape=21,
             aes(x=Dim.1, y=Dim.2, 
                 fill = city,
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
  xlab(paste("PC1 (", round(D1VE, 2), "% VE)" , sep = "")) +
  ylab(paste("PC2 (", round(D2VE, 2), "% VE)", sep = "")) ->
  PCA_12 


#Save the plot
ggsave(PCA_12, 
       file = paste("PCA", "pdf" , sep = ".") ,
       width = 4,
       height = 4)

#############
lm(Dim.2 ~ collectionDate + lat + long + ,
   data = PCA_table) %>% summary()



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