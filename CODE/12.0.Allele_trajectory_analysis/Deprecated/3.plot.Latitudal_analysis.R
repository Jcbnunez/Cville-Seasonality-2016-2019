### libraries

library(tidyverse)
library(magrittr)
library(data.table)
library(car)
library(DescTools)
library(foreach)
library(doMC)
library(patchwork)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(lubridate)
library(forcats)
library(viridis)
library(SeqArray)
library(tidyverse)
library(vcfR)
library(scatterpie)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)

registerDoMC(2)

setwd("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/")


### load meta-data file
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

### common SNP.dt
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))
snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, snp.dt$id)

### function
getData <- function(chr="2L", start=14617051, end=14617051) {
  # chr="2L"; start=14617051; end=14617051
  
  ### filter to target
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)
  
  ### get annotations
  #message("Annotations")
  tmp <- seqGetData(genofile, "annotation/info/ANN")
  len1 <- tmp$length
  len2 <- tmp$data
  
  snp.dt1 <- data.table(len=rep(len1, times=len1),
                        ann=len2,
                        id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))
  
  # Extract data between the 2nd and third | symbol
  snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
  snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
  
  # Collapse additional annotations to original SNP vector length
  snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                        list(variant.id=id)]
  
  snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
  snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  afi <- merge(af, snp.dt1.an, by="variant.id")
  afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis[,-c("n"), with=F]
}

### test
data <- getData()



#load latitude model
load("./latitude_model.Rdata")

#load temperature model
glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"
load(glm.file)
glm.out %>%
  filter(mod == "aveTemp+year_factor",
         chr == "2L" ) %>% 
  mutate(window = case_when(
    pos > 5155762 & pos < 5255762 ~ "win5",
    pos > 6255762 & pos < 6355762 ~ "win6",
    pos > 9505762 & pos < 9605762 ~ "win9")) %>%
  filter(window %in% c("win5", "win6", "win9")) -> outliers


load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"

padding=0
starts= c(5155762, 6255762, 9505762)
ends= c(5255762, 6355762, 9605762)

win_5np = getData(chr="2L", start=5155762-padding , end=5255762+padding) 
win_6np = getData(chr="2L", start=6255762-padding , end=6355762+padding) 
win_9np = getData(chr="2L", start=9505762-padding , end=9605762+padding)

### join datasets
rbind(mutate(win_5np, win = "1.win5"), 
      mutate(win_6np, win = "2.win6"),
      mutate(win_9np, win = "3.win9")) %>%   
  left_join(weather.ave) ->
  win_outliers_no_padding

### Make additional files;, simulans polarization and annotations
win_outliers_no_padding %>%
  filter( set %in% c("dgn"),
          sampleId == "SIM") %>%
  dplyr::select(chr, pos, sim_af = af) -> sim_af

win_outliers_no_padding %>%
  group_by(pos) %>%
  slice_head(n=1) %>%
  dplyr::select(pos, col) -> annots

####
outliers %<>%
  separate(b, into = c("b1","b2","b3","b4","b_temperature"), sep = ";") 

left_join(outliers, latitude_model) %>% 
 .[complete.cases(.$mod),] %>% 
  full_join(annots) %>% 
  mutate(window = RoundTo(pos, 1e6, "floor"),
         concordance = case_when(sign(as.numeric(b_temperature)) != sign(lat_beta) ~ 1,
                                 sign(as.numeric(b_temperature)) == sign(lat_beta) ~ -1),
         func_priority = case_when(col %in% 
                                     c("missense_variant", 
                                       "3_prime_UTR_variant", 
                                       "5_prime_UTR_variant") ~ "yes",
                                   !(col %in% 
                                       c("missense_variant", 
                                         "3_prime_UTR_variant", 
                                         "5_prime_UTR_variant")) ~ "no")
                                         ) ->
  lat_seasonal_concordance

lat_seasonal_concordance %>%
  filter(
        #!is.na(concordance),
         p.lrt < 0.05) %>% 
     ggplot(
    aes(
      x=pos/1e6,
      y= -log10(p.adjust(lat_p_val, "fdr"))*concordance,
      color=as.factor(concordance),
      fill=as.factor(concordance),
      #shape = func_priority
    )
  ) + geom_point(size = 1)+
  geom_hline(yintercept = -log10(0.05)) +
  geom_hline(yintercept = log10(0.05)) +
  scale_shape_manual(values = c(21,22)) +
  ylab("Lat cor P * (Temp sign)") +
  xlab("Genomic Position (Mb)") +
  facet_wrap(~window, scales = "free_x") ->
  betas_plot

ggsave(betas_plot, file = "betas_plot.pdf", w= 8, h=2.3)

####
####
### 
### 

latitude_model %>%
  mutate(window = RoundTo(pos, 1e6, "floor") ) %>%
  ggplot(aes(
    x=pos,
    y=-log10(lat_p_val),
    color =lat_beta
  )) +
  geom_vline(data = data.frame(pos = c(5160331, 6354995, 9576345), window = c(5e6, 6e6, 9e6)), 
             aes(xintercept = pos)) +
  geom_point() +
  scale_color_gradient2() +
  facet_wrap(~window, scales = "free_x", ncol = 1) ->
  lat_model_plot
ggsave(lat_model_plot, file = "lat_model_plot.pdf")
         
### Plot in a world map
win_outliers_no_padding_DEST_polarized %>% 
  .[complete.cases(.$af_neff_pol),] %>%
  filter(pos == 5160331,
         set %in% c(#"DrosRTEC",
                    "DrosEU"),
         year > 2012,
         lat > 30 & lat < 66,
         ) %>%
  group_by(locality, lat, long, country) %>%
  summarize(Ancestral = mean(af_neff_pol),
            Derived = 1-mean(af_neff_pol)) -> Msp300_EU

## create world maps
#Graph maps
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(
    #xlim = c(-125.15, 40.00), ylim = c(20, 65.00), 
    xlim =  c(-12, 41.00), ylim = c(38.00, 63.00),
    #xlim =  c(-99, -43), ylim = c(5, 65),
    expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) -> Europe

Europe + 
  geom_scatterpie(data =Msp300_EU , 
             aes(x=long, 
                 y=lat, 
                 group=locality,
                 r = 1.1
                 ),
             cols=c("Ancestral", "Derived"),
             size = 0.6,  alpha = 0.7, #shape = 21
             )  +
  xlab("Lon") + ylab("Lat") + ggtitle("Clin. of Msp300-5160331 ") + 
  #facet_wrap(~year, ncol= 1) +
  theme(legend.position = "top") ->
  AF_pie_all

ggplot(data = Msp300_EU,
       aes(x=as.factor(RoundTo(lat, 5 , "floor")),
           y=Derived)) +
  geom_boxplot() +
  #facet_wrap(~year, ncol= 1) +
  coord_flip() +
  ylab("Derived AF") + xlab("Frequency") +
  theme(legend.position = "top") ->
  AF_box_all
###
EW_dat = filter(Msp300_EU, !(country %in% c("United Kingdom", "Serbia","Ukraine", "Turkey", "Hungary" ) ))
       
Europe + 
  geom_scatterpie(data = EW_dat ,
                  aes(x=long, 
                      y=lat, 
                      group=locality,
                      r = 1.1),
                  cols=c("Ancestral", "Derived"),
                  size = 0.6,  alpha = 0.7)  +
  xlab("Lon") + ylab("Lat") + ggtitle("Clin. of Msp300-5160331 ") + 
  #facet_wrap(~year, ncol= 1) +
  theme(legend.position = "top") ->
  AF_pie_cline_EW

ggplot(data = EW_dat ,
       aes(x=as.factor(RoundTo(lat, 5 , "floor")),
           y=Derived)) +
  geom_boxplot() +
  #facet_wrap(~year, ncol= 1) +
  coord_flip() +
  ylab("Derived AF") + xlab("Frequency") +
  theme(legend.position = "top") ->
  AF_box_cline_EW
#####
###
NS_dat = filter(Msp300_EU, country %in% c("United Kingdom", "France") )

Europe + 
  geom_scatterpie(data = NS_dat ,
                  aes(x=long, 
                      y=lat, 
                      group=locality,
                      r = 1.1),
                  cols=c("Ancestral", "Derived"),
                  size = 0.6,  alpha = 0.7)  +
  xlab("Lon") + ylab("Lat") + ggtitle("Clin. of Msp300-5160331 ") + 
  #facet_wrap(~year, ncol= 1) +
  theme(legend.position = "top") ->
  AF_pie_cline_NS

ggplot(data = NS_dat ,
       aes(x=as.factor(RoundTo(lat, 5 , "floor")),
           y=Derived)) +
  geom_boxplot() +
  #facet_wrap(~year, ncol= 1) +
  coord_flip() +
  ylab("Derived AF") + xlab("Frequency") +
  theme(legend.position = "top") ->
  AF_box_cline_NS

####
ggsave((AF_pie_all+AF_box_all)/
       (AF_pie_cline_EW+AF_box_cline_EW)/ 
       (AF_pie_cline_NS+AF_box_cline_NS),
       file = "clinality.pdf", w=6, h=9)



##### NORTH AMERRICA
##### NORTH AMERRICA
##### NORTH AMERRICA
##### NORTH AMERRICA
##### NORTH AMERRICA
##### NORTH AMERRICA
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(
    #xlim = c(-125.15, 40.00), ylim = c(20, 65.00), 
    #xlim =  c(-12, 41.00), ylim = c(38.00, 63.00),
    xlim =  c(-126, -70), ylim = c(30, 65),
    expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) -> AmericaN


win_outliers_no_padding_DEST_polarized %>% 
  .[complete.cases(.$af_neff_pol),] %>%
  filter(pos == 5160331,
         continent %in% c("North_America", "NorthAmerica"),
         year > 2012,
  ) %>%
  group_by(locality, lat, long, country) %>%
  summarize(Ancestral = mean(af_neff_pol),
            Derived = 1-mean(af_neff_pol)) -> Msp300_AmN

East = filter(Msp300_AmN, long > -100  )
West = filter(Msp300_AmN, long < -100  )

AmericaN + 
  geom_scatterpie(data = East  , 
                  aes(x=long, 
                      y=lat, 
                      group=locality,
                      r = 1.1
                  ),
                  cols=c("Ancestral", "Derived"),
                  size = 0.6,  alpha = 0.7, #shape = 21
  )  +
  xlab("Lon") + ylab("Lat") + ggtitle("Clin. of Msp300-5160331 ") + 
  #facet_wrap(~year, ncol= 1) +
  theme(legend.position = "top") ->
  AF_pie_AmericaE

ggplot(data = East,
       aes(x=as.factor(RoundTo(lat, 5 , "floor")),
           y=Derived)) +
  geom_boxplot() +
  #facet_wrap(~year, ncol= 1) +
  coord_flip() +
  ylab("Derived AF") + xlab("Lat") +
  theme(legend.position = "top") ->
  AF_box_AmericaE



AmericaN + 
  geom_scatterpie(data = West  , 
                  aes(x=long, 
                      y=lat, 
                      group=locality,
                      r = 1.1
                  ),
                  cols=c("Ancestral", "Derived"),
                  size = 0.6,  alpha = 0.7, #shape = 21
  )  +
  xlab("Lon") + ylab("Lat") + ggtitle("Clin. of Msp300-5160331 ") + 
  #facet_wrap(~year, ncol= 1) +
  theme(legend.position = "top") ->
  AF_pie_AmericaW

ggplot(data = West,
       aes(x=as.factor(RoundTo(lat, 5, "floor")),
           y=Derived)) +
  geom_boxplot() +
  #facet_wrap(~year, ncol= 1) +
  coord_flip() +
  ylab("Derived AF") + xlab("Lat") +
  theme(legend.position = "top") ->
  AF_box_AmericaW
####
ggsave((AF_pie_AmericaE+AF_box_AmericaE)/
        (AF_pie_AmericaW+AF_box_AmericaW),
       file = "clinality.NAM.pdf", w=8, h=8)

