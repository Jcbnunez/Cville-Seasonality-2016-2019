
rm(list = ls())
# Load packages

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
library(devtools)
library(lubridate)


load("make.Fig2.panelA.Rdata")

ggplot() + 
  geom_point(shape =21, 
             size = 3,
             alpha = 0.9,
             data = pca_table_df,
             aes(
               x=Dim.1,
               y=Dim.2,
               fill = year)
  ) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 2016) +
  theme_bw() +
  ylab("PC2") +
  xlab("PC1") +
  facet_wrap(~chr) ->
  cville_chr_pca

ggsave(cville_chr_pca,
       file = "cville_chr_pca.pdf",
       width = 5,
       height = 4)

### Quantify
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "2L" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "2R" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "3L" )))
anova(lm(Dim.1 ~ year, data = filter(pca_table_df, chr == "3R" )))


anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "2L" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "2R" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "3L" )))
anova(lm(Dim.1 ~ `In(2L)t`, data = filter(pca_table_df, chr == "3R" )))


