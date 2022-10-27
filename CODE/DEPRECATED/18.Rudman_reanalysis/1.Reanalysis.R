####
####
####
####
library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(gmodels)
library(doMC)
registerDoMC(4)
library(tidyr)

files.to.load = 
system("ls /project/berglandlab/Dmel_genomic_resources/Datasets/Experimental_Evolution_Rudman2022/ | grep '2L' ",
       intern = T)

rudman.afs = 
foreach(i=1:length(files.to.load), .combine = "rbind")%do%{
  
  message(paste(i,"/",length(files.to.load)))
  
  metadata = tstrsplit(files.to.load[i], "\\.")
  
  fread(paste("/project/berglandlab/Dmel_genomic_resources/Datasets/Experimental_Evolution_Rudman2022/",
              files.to.load[i], 
              sep = "" )) %>%
    mutate(name = metadata[[1]],
           chr = metadata[[5]]) 
}

rudman.afs %<>%
  separate(name, into = c("Cage", "Sex", "Timepoint"), sep = "_" )

save(rudman.afs, file = "rudman.raw.afs.Rdata")

rudman.afs %<>%
  mutate(win = 
           case_when(
             pos > 4650065 & pos < 4799922 ~ "win_4.7",
             pos > 5100324 & pos < 5349218 ~ "win_5.2",
             pos > 6100321 & pos < 6349489 ~ "win_6.2",
             pos > 9500286 & pos < 9700005 ~ "win_9.6"
           ))
###
load("haplo_tags_SNPids.Rdata")
haplo_tags_SNPids %>%
  separate(SNP_id, into = c("chr", "pos", "type"), sep = "_" ) %>%
  mutate(pos = as.numeric(pos))  ->
  haplotags

###
rudman.afs %>%
  filter(Timepoint > 1) %>%
  filter(!is.na(win)) %>%
  filter(af > 0.01) %>%
  mutate(pos = as.numeric(pos)) %>%
  right_join(haplotags) %>%
  filter(!is.na(win)) %>%
  group_by(Cage, Sex, Timepoint, chr,  win) %>%
  summarize(af_me = mean(af),
            af_lci = ci(af)[2],
            af_uci = ci(af)[3]
            ) %>%
  ggplot(
    aes(
      x=as.numeric(Timepoint),
      y=as.numeric(af_me),
      ymax=af_lci,
      ymin=af_uci,
      #color = as.factor(pos)
    )
  ) + 
  geom_errorbar(width = 0.1) +
  geom_point(aes(color = Cage)) +
  geom_smooth(color = "black") +
  #theme(legend.pos = "none") +
  facet_grid(.~win) ->
  rudman_plot

ggsave(rudman_plot, file = "rudman_plot.png", w = 5, h = 3)


