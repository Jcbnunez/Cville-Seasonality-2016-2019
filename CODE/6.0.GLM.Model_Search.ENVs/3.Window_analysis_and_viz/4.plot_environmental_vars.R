### weather
### 
### 
### 

library(tidyverse)
library(data.table)
library(patchwork)

demo.assg <- fread("/project/berglandlab/Dmel_genomic_resources/Datasets/DEST.v.1.Demography/DEST_Sample_clusters.txt")
samp.data <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")


load("/project/berglandlab/Dmel_genomic_resources/Datasets/NASA_power_weather/nasa_power.weather.mod.Rdata")

dim(weather.ave)[2] - 2
weather.ave$mod %>% table  
10*11

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  full_join(demo.assg) %>%
  full_join(dplyr::select(samp.data,sampleId, collectionDate, country, year, locality, collectionDate) ) %>%
  mutate(as.date.col = as.Date(collectionDate, format = "%m/%d/%Y"))->
  weather.ave.meta


weather.ave.meta %>%
  filter(Continental_clusters == "3.Europe_E" & mod == 9) %>%
  ggplot(
    aes(
      x=as.date.col,
      y=temp.ave,
      color = as.factor(year)
    )
  ) + geom_point() +
  geom_smooth() +
  ggtitle("Temp 9 EU E") +
  facet_grid(.~year, scales = "free_x")->
  t.ave.9.EU_E

weather.ave.meta %>%
  filter(Continental_clusters == "1.Europe_W" & mod == 8) %>%
  ggplot(
    aes(
      x=as.date.col,
      y=humidity.ave,
      color = as.factor(year)
    )
  ) + geom_point() +
  geom_smooth() +
  ggtitle("Temp 9 EU W") +
  facet_grid(.~year, scales = "free_x")->
  t.humidity.8.EU_W


weather.ave.meta %>%
.[grep("VA_ch", .$sampleId),] %>%
filter( mod == 2) %>%
  separate(sampleId, into = c("State", "City", "Year", "Date"), sep = "_" , remove = F) %>%
filter(Year > 15) %>%
  mutate(complex.dat = paste(com.dat = paste(Date, Year, sep = "y"))) %>%
  mutate(as.date.col = as.Date(complex.dat, format = "m%md%dy%Y")) %>%
  ggplot(
    aes(
      x=as.date.col,
      y=temp.max,
      color = as.factor(Year)
    )
  ) + geom_point() +
  geom_smooth() +
  ggtitle("Temp Max 2 Cville") +
  facet_grid(.~Year, scales = "free_x")->
  t.max.2.cville

ggsave(t.ave.9.EU_E/ t.humidity.8.EU_W / t.max.2.cville, file = "cville.europe.best.models.pdf", h = 6, w = 9)

