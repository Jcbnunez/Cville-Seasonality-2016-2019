## Process downloaded material from Flybase
## make sure to request a batch download with the following data:
## SUBMITTED ID	GENE SNAPSHOT	LOCATION ARM	LOCATION MAX	LOCATION MIN	LOCATION STRAND	NAME	SYMBOL
## 

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)


# import the flybase portion
FlyBase_Fields_download <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/8.Gowinda/human_curated/FlyBase_Fields_download.txt")
FlyBase_Fields_download %<>%
  mutate(Median_pos = (as.numeric(LOCATION_MAX)+as.numeric(LOCATION_MIN))/2)


FlyBase_Fields_download %>% 
  .[complete.cases(.$Median_pos),] %>% 
  group_by(SYMBOL) %>%
  slice_head() ->
  FlyBase_Fields_download_dedup


## Now import the go portion
GOwinda_curated_of_interest_bimin5 <- read.delim("~/Documents/GitHub/Cville-Seasonality-2016-2019/8.Gowinda/human_curated/GOwinda_curated_of_interest_bimin5.txt", sep = "\t")
gene_list=list()
for(i in 1:dim(GOwinda_curated_of_interest_bimin5)[1]){
  
  
  gene_list[[i]] = data.frame(X.SUBMITTED.ID = 
               toupper(unlist(strsplit(GOwinda_curated_of_interest_bimin5$Genes[i], split=",")[[1]])),
             FDR= GOwinda_curated_of_interest_bimin5$FDR[i],
             Nobs= GOwinda_curated_of_interest_bimin5$Nobs[i],
             GO.id = GOwinda_curated_of_interest_bimin5$GO.id[i],
             GO.term= GOwinda_curated_of_interest_bimin5$GO.term[i],
             Supercategory = GOwinda_curated_of_interest_bimin5$Supercategory[i],
             Supracategory = GOwinda_curated_of_interest_bimin5$Supracategory[i]
             )
    
}

gene_df = do.call(rbind, gene_list)

left_join(FlyBase_Fields_download_dedup, gene_df ) ->
  GO_fly_base_dat
  
GO_fly_base_dat %>%
  ggplot(aes(
    x=Median_pos,
    y=GO.id,
             )) +
  geom_point(alpha = 0.1) +
  facet_wrap(~Supercategory, scale = "free")  +
  geom_vline(xintercept = 2225744) + 
  geom_vline(xintercept = 13154180) 


#### Count genes per window
#### Make windows
win.bp <- 1e5
step.bp <- 5e4

setDT(GO_fly_base_dat)
setkey(GO_fly_base_dat, "Median_pos")

tmp <- GO_fly_base_dat
wins = data.table(
  start=seq(from=min(tmp$Median_pos), to=max(tmp$Median_pos)-win.bp, by=step.bp),
  end=seq(from=min(tmp$Median_pos), to=max(tmp$Median_pos)-win.bp, by=step.bp) + win.bp)

wins[,i:=1:dim(wins)[1]]
dim(wins)

setkey(GO_fly_base_dat, "Median_pos")

  win.out <- foreach(win.i=c(1:dim(wins)[1]), 
                     .combine = "rbind",
                     .errorhandling="remove")%dopar%{
                       
                       win.tmp <- GO_fly_base_dat %>%
                         filter(Median_pos >= wins$start[win.i],
                                Median_pos <= wins$end[win.i],
                         )
                       
                       win.tmp %>%
                         as.data.frame() %>%
                         group_by(GO.id) %>%
                         summarise(
                           window=median(wins$start[win.i],wins$end[win.i]),
                           N = n())
                       
                     }## close do par

  win.out %>%
    ggplot(
      aes(
        x=window,
        y=GO.id,
        fill=N,
      ) 
    ) +
    geom_raster() +
    geom_vline(xintercept = 2225744) + 
    geom_vline(xintercept = 13154180) +
    theme_classic() +
    ggtitle("Number of Genes at a window | conditioned by GO category") +
    scale_fill_gradient2(low = "blue", mid = "steelblue", high = "red", midpoint = 3)
  

#### Case Studies

  GO_fly_base_dat %>%
    filter(Median_pos > 10330750-1e5,
           Median_pos < 10330750+1e5,
           GO.id %in% c("GO:0044422"))


  
  
