### sample selection for spatial analysis
### 
library(data.table)
library(tidyverse)
library(foreach)
library(SeqArray)
library(gmodels)

library(rnaturalearth)
library(rnaturalearthdata)

######
## We will use generated in step 4
load("../12.trajectory_analysis/haplo_tags_SNPids.Rdata")
haplo_tags_SNPids

haplo_tags_SNPids %<>% 
  filter(class == "GLM_LD") %>%
  separate(SNP_id, into = c("chr","pos", "type"), remove = F) %>%
  mutate(win = 
           case_when(
             pos > 4650065 & pos < 4799922 ~ "win_4.7",
             pos > 5100324 & pos < 5349218 ~ "win_5.2",
             pos > 6100321 & pos < 6349489 ~ "win_6.2",
             pos > 9500286 & pos < 9700005 ~ "win_9.6"
           ))  
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
final_in2Lt_markers %<>%
  data.frame(SNP_id = .) %>%
  separate(SNP_id, into = c("chr", "pos", "type"), sep = "_", remove = F) %>%
  mutate(win = "inv",
         class = "inv")

rbind(haplo_tags_SNPids, final_in2Lt_markers) -> haplo_tags_SNPids_and_inv





#######
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

### get haplotags
haplotag_snps_AFS = foreach(i=1:dim(haplo_tags_SNPids_and_inv)[1], .combine = "rbind")%do%{
  
  snp_tmp <- getData(chr="2L", 
                     start=haplo_tags_SNPids_and_inv$pos[i], 
                     end=haplo_tags_SNPids_and_inv$pos[i]) %>%
    mutate(win = haplo_tags_SNPids_and_inv$win[i]  )
  return(snp_tmp)
}

#### Polarize
haplotag_snps_AFS %>%
  filter(sampleId == "SIM") %>%
  dplyr::select(variant.id, pos, sim_af=af) ->
  SIM_AF

left_join(haplotag_snps_AFS, SIM_AF) %>%
  mutate(af_polarized = case_when(
    sim_af == 0 ~ 1-as.numeric(af),
    sim_af == 1 ~ as.numeric(af))) %>%
  mutate(collectionDate = as.Date(collectionDate, format = "%m/%d/%Y")) ->
  haplotag_snps_AFS_pol
####
####
####
####
####
####




haplotag_snps_AFS_pol %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC") ) %>%
  filter(!is.na(af_polarized)) %>%
  group_by(sampleId, collectionDate, set, year, win, locality, continent, season, lat, long) %>%
  summarise(Mean_haplotag = mean(af, na.rm = T),
            sd_haplotag = sd(af, na.rm = T),
            ci_l = ci(af, na.rm = T)[2],
            ci_h = ci(af, na.rm = T)[3],
            locality = locality,
            year = year,
            season = season ,
            continent = continent,
            lat = lat,
            long = long
  ) -> haplotag_snps_AFS_sums
  

haplotag_snps_AFS_sums %>%
  filter(locality != "VA_ch") %>%
  filter(continent == "NorthAmerica") %>%
  filter(year %in% 2014) %>% 
  filter(season == "fall") -> nam_pts

###
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

###
ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-100.15, -50.0), 
           ylim = c(23.00, 52.00), 
           expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(data =nam_pts , 
             aes(x=as.numeric(long), 
                 y=as.numeric(lat), 
                 fill = scale(Mean_haplotag), 
                 #shape = as.factor(Continental_clusters)
                 ),
             size = 2.5,  
             alpha = 0.7,
             shape  = 21)  +
  scale_fill_gradient2(low = "steelblue", high = "firebrick", na.value = NA, 
                       #midpoint = 0.5
                       ) + 
  ggtitle("2014 Fall samps NAmerica") +
  xlab("Lon") +
  ylab("Lat") + 
  #theme(legend.position = "none") + 
  facet_wrap(~win) -> NAmerica_pts

ggsave(NAmerica_pts, file = "NAmerica_pts.pdf")


###

haplotag_snps_AFS_sums %>%
  filter(locality != "VA_ch") %>%
  filter(continent == "Europe") %>%
  filter(year == 2014) %>% 
  filter(season == "fall") -> Eu_pts


ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim =  c(-12, 41.00), ylim = c(32.00, 63.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_point(data =Eu_pts , 
             aes(x=as.numeric(long), 
                 y=as.numeric(lat), 
                 fill = scale(Mean_haplotag), 
                 #shape = as.factor(Continental_clusters)
             ),
             size = 2.5,  
             alpha = 0.7,
             shape  = 21)  +
  scale_fill_gradient2(low = "steelblue", high = "firebrick", na.value = NA, 
                       #midpoint = 0.5
  ) + 
  ggtitle("2014 Fall samps EU") +
  xlab("Lon") +
  ylab("Lat") +
  facet_wrap(~win) -> 
  Eu_pts_fig

ggsave(Eu_pts_fig, file = "Eu_pts_fig.pdf")

#####
haplotag_snps_AFS_sums %>%
  filter(locality != "VA_ch") %>%
  #filter(continent == "NorthAmerica") %>%
  filter(year %in% 2014) %>% 
  filter(season == "fall") %>% 
  ggplot(
    aes(
      x=lat,
      y=Mean_haplotag
    )
  ) + geom_point() +
  facet_grid(win~continent, scales = "free_x") ->
  lat_plot

ggsave(lat_plot, file = "lat_plot.pdf")


haplotag_snps_AFS_sums %>%
  filter(locality != "VA_ch") %>%
  #filter(continent == "NorthAmerica") %>%
  filter(year %in% 2014) %>% 
  filter(season == "fall") %>% 
  ggplot(
    aes(
      x=long,
      y=Mean_haplotag
    )
  ) + geom_point() +
  facet_grid(win~continent, scales = "free_x") ->
  long_plot

ggsave(long_plot, file = "long_plot.pdf")

