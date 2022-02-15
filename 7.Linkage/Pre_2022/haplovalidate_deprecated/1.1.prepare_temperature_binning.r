#### Prepare the data 
#### 

library(data.table)
library(tidyverse)
library(reshape2)
library(magrittr)
library(gmodels)

##load in the data
load("/project/berglandlab/jcbnunez/Shared_w_Connor/Drosophila/Cville_2L.ECfiltered.Rdata")
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")

names(weather.ave)[1] = "sampleId"

filtered_samps_for_analysis %>%
  filter(MeanEC > 30,
         set == "CvilleSet") %>% 
  left_join(weather.ave) -> data_for_slice

####
####
#data_for_slice %>%
#  filter(Month == 10)

library(factoextra)
data_for_slice %>%
  select(yday,
         aveTemp) -> dat_in

  km.res <- kmeans(scale(dat_in), 10, nstart = 25)

  
  cbind(data_for_slice, cluster = km.res$cluster) %>%
    ggplot(
      aes(
        x=yday,
        y=aveTemp,
        color = as.factor(cluster),
        shape = as.factor(year),
        label = sampleId
      )
    ) +
    geom_text(size = 2) +
    geom_point() -> gg_cluster
  
    ggsave(gg_cluster,
           file = "gg_cluster.pdf")
  
    
    clusters <- c(
      
      "VA_ch_18_m07d12", #State 1
      "VA_ch_17_m07d07", #State 1
      "VA_ch_16_m07d08", #State 1
    
      "VA_ch_16_m09d02", #State 2
      "VA_ch_17_m08d17", #State 2
      "VA_ch_18_m09d06", #State 3
      
      "VA_ch_16_m09d16", #State 4
      "VA_ch_18_m09d20", #State 4
      "VA_ch_17_m10d12", #State 4
      
      "VA_ch_18_m10d18", #State 5
      "VA_ch_16_m10d14", #State 5
      "VA_ch_17_m10d26", #State 5
      
      "VA_ch_16_m11d11", #State 6
      "VA_ch_18_m11d01", #State 6
      "VA_ch_17_m11d09"  #State 6
    )
    
data_for_slice %>%
  filter(sampleId %in%
           clusters
  ) %>%
  ggplot(aes(
    x=aveTemp,
    y=year,
    color = as.factor(year)
  )) + 
  geom_point() ->
  plot_year_v_temp

ggsave(plot_year_v_temp,
       file = "plot_year_v_temp.pdf")
