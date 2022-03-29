library(directlabels)
library(googlesheets4)
library(data.table)
library(censusxy)
library(foreach)
library(elevatr)
library(ggmap)
library(ggplot2)
library(viridis)
library(cowplot)
library(elevatr)
library(rnoaa)
library(maps)
library(data.table)
library(raster)
library(patchwork)
library(forcats)
library(ggtree)
library(ape)
library(reshape2)
library(tidyverse)

#Aim 1: phylogenetic prep and species abundance.
#

Orchard_counts <- read.delim("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/2.1.Species_count_data/speciesCounts.txt")

Orchard_counts$species = gsub("Chymomyza amoena", "other", Orchard_counts$species)
Orchard_counts$species = gsub("unknown", "other", Orchard_counts$species)
Orchard_counts$species = gsub("Drosophila tripunctata", "other", Orchard_counts$species)
Orchard_counts$species = gsub("Drosophila mel/sim", "other", Orchard_counts$species)

Orchard_counts$species = gsub("Drosophila", "", Orchard_counts$species)
Orchard_counts$species = gsub(" indianus", "", Orchard_counts$species)


colors = c("red",	
           "blue",	
           "green",	
           "yellow",	
           "purple",	
           "orange",
           "grey")

Orchard_counts %>% 
  separate(CollectionDate,
           into = c("Year", "MD"), sep = 4) %>%
  separate(MD,
           into = c("Month", "Day"), sep = 2) %>%
  group_by(Month, species, Year) %>%  
  summarise(Species_count = sum(Count))  %>%
  mutate(habitat = "Orchard") ->
  orchatrd_dat

orchatrd_dat$species = gsub("^ ", "", orchatrd_dat$species)
orchatrd_dat$species = as.factor(orchatrd_dat$species)
orchatrd_dat$Month = as.integer(orchatrd_dat$Month )
orchatrd_dat$Year = as.double(orchatrd_dat$Year )

### Compost pile
Collections_2020season <- read.delim("~/Dropbox/2021_BackEvo/Collections_2020season.txt")
Collections_2020season$Date = as.Date(Collections_2020season$Date, format = "%m/%d/%Y" )
Collections_2020season$Collector = gsub("KrSa","FrSa", Collections_2020season$Collector )

Collections_2020season %>%
  .[,which(names(.) %in% c(  "melanogaster",	
                             "suzukii",	
                             "immigrans",	
                             "simulans",	
                             "hydei",	
                             "zaprionus", 
                             "other",
                             "Collector", 
                             "Tube.Number", 
                             "Date",  
                             "Date.Sorted", 
                             "Sorted.by"))] %>%
  melt(id =c("Collector", 
             "Tube.Number", 
             "Date",  
             "Date.Sorted", 
             "Sorted.by"))  %>%
  mutate(julian_date = yday(Date),
         Month = month(Date)) -> Collections_2020season_m

Collections_2020season_m %>%
  .[which(.$Month != 12),] %>%
  group_by(variable, Month) %>% 
  summarise(Species_count = sum(value)) %>%
  mutate(habitat = "Compost",
         Year = 2020) ->
  compost_pile_dat
names(compost_pile_dat)[1] = "species"
compost_pile_dat$species = gsub("zaprionus", "Zaprionus", compost_pile_dat$species)

rbind(orchatrd_dat, compost_pile_dat) %>%
  group_by(species, habitat) %>%
  summarise(Nflies = sum(Species_count)) %>%
  dcast(species~habitat, value.var = "Nflies")


colors2 = c("red",	
           "blue",	
           "green",	
           "yellow",	
           "purple",	
           "orange")


rbind(orchatrd_dat, compost_pile_dat) %>%
  .[which(.$species != "other"),] %>%
  group_by(species, Year, Month, habitat) %>%
  summarise(Nflies = sum(Species_count)) %>% 
  ggplot(aes(
    x=Month,
    y=Nflies,
    color = species,
    linetype = as.factor(Year),
    shape = as.factor(Year)
  )) +
  geom_smooth(se = F,
              span = 1,
              size = 0.5) +
  geom_point(size = 2) +
  facet_wrap(habitat~species, 
             scales = "free",
             nrow = 2) +
  theme_dark()  +
  scale_color_manual(values = colors2) +
  theme(legend.position = "top" )

### Only Melanogaster  

rbind(orchatrd_dat, compost_pile_dat) %>%
  filter(species == "melanogaster",
         Year != 2020) %>%
  #left_join(N_tot_flies) %>%
  ggplot(aes(
    x=Month,
    y=Species_count,
    color = as.factor(Year),
    linetype = as.factor(Year),
    shape = as.factor(Year)
  )) +
  geom_smooth(se = F,
              span = 10,
              size = 0.5) +
  geom_point(size = 3.5) +
  ylab("Indviduals") +
  theme_bw()  +
  theme(legend.position = "top" )






################  
rbind(orchatrd_dat, compost_pile_dat) %>%
 ggplot(aes(x=as.factor(Month),
           y=Species_count,
           fill = gsub("other", "zu.other" ,species)
)) +
  geom_bar(stat = "identity",
           position="fill",
           alpha = 1)  +
  theme_cowplot() + 
  scale_fill_manual(values = colors) +
  facet_wrap(Year~habitat, scales = "free_x") + 
  labs(fill = "Species") + 
  theme(legend.position="top") +
  ylab("Collection Month") +
  xlab("Abundance (%)")



####
####
# Make Tress
####
####


traits <- read.delim("~/Documents/JCBN_files/Grants/2021_AOB_career/Data/Aim1_phylo/traits.txt")
nwk <- "/Users/jcbnunez/Documents/JCBN_files/Grants/2021_AOB_career/Data/Aim1_phylo/drosophilid_fasta.aln.fasta.contree"
tree <- read.tree(nwk)
nodelabels()

nwk_tr <- root(tree, 7)
plot(nwk_tr)

tree <- full_join(nwk_tr, traits, by = 'tip.label')

ggtree(nwk_tr,
       layout="slanted",
       branch.length = 'none',
       size = 2)  %<+% 
  traits  ->
  tree_obj_traits


tree_obj_traits + 
  geom_tiplab() + 
  geom_rootedge() +
  xlim(-1, 9)


tree_obj_traits$data %>%
  as.data.frame() %>%
  .[which(.$isTip == T),] %>%
  ggplot(aes(
    x=x,
    ymin=y-0.1,
    ymax=y+0.1,
    xmin=CTmin,
    xmax=CTmax,
  )) +
  geom_rect( aes(fill=),fill="red") +
  geom_vline(xintercept = 25, linetype = "dashed") +
  theme_classic()


#####
#####
#####

