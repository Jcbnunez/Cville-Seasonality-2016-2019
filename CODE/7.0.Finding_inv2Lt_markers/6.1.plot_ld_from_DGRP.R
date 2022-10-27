# Plot the LD analysis from the DGRP

library(data.table)
library(tidyverse)
library(magrittr)
library(reshape2)
library(gmodels)

#define 2lt region
in2lt_beg=2225744	
in2lt_end=13154180

#load in data
load("./ld_df.rank90_99.Rdata")

#modify rank variable
ld_df$rank = as.numeric(ld_df$rank)

## summarize
ld_df %>%
  group_by(rank) %>%
  summarise(meanR2 = mean(R2),
            sdR2 = sd(R2)
            )

ld_df %>%
  filter(rank == 0.99, R2 > 0.7) %>% 
  ggplot(
    aes(
    x=as.numeric(BP_B),
    y=as.numeric(R2),
    color = as.numeric(rank)
  )
  ) + 
  geom_point() +
  geom_vline(xintercept = in2lt_beg) +
  geom_vline(xintercept = in2lt_end) +
  facet_wrap(~id_tag)->
  r2_dgrp

ggsave(r2_dgrp,
       file = "r2_dgrp.png",
       width = 20,
       height = 20
       )
  
  
  
  