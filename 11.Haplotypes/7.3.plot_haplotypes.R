library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)

load("./snp_dat_table_loc_homkar_melt_pos_samps.Rdata")

snp_dat_table_loc_homkar_melt_pos_samps$pos = as.numeric(snp_dat_table_loc_homkar_melt_pos_samps$pos)

snp_dat_table_loc_homkar_melt_pos_samps %>%
  ggplot(
    aes(
      x=as.numeric(pos_id),
      y=samp_hap,
      fill = as.factor(value)
    )
  ) + geom_tile() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  facet_wrap(~karyo, ncol = 1, scales = "free") ->
  karyo_plot

ggsave(karyo_plot, file = "karyo_plot.png")
