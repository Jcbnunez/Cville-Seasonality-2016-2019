### Genetics and Environemtns characterization of Genetic targets

## Part 1 -- Identify inversion status

library(tidyverse)
library(magrittr)
library(data.table)
library(car)
library(DescTools)

## Load SVM predictions
## uses ../Cville-Seasonality-2016-2019/5.Finding_inv2Lt_markers/9.train_predictive_model.r

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")

CM_predictions %>%
  ggplot(aes(SVM)) +
  geom_vline(xintercept = 0.90) +
  geom_vline(xintercept = 0.10) +
  geom_histogram() -> svm_hist

ggsave(svm_hist, file = "svm_hist.pdf")

CM_predictions %>%
  separate(line_name, into = c("set", "id", "date"), sep = "_") %>%
  mutate(date_parsed = case_when( date %in% c("CMspring") ~ as.Date("06/01/2018", format = "%m/%d/%Y"),
                                  date %in% c("CMfall") ~ as.Date("09/01/2017", format = "%m/%d/%Y"),
                                  !(date %in% c("CMspring", "CMfall")) ~ as.Date(paste(date, 2016, sep = ""), format = "%m%d%Y"))) %>%
  mutate(karyot = case_when(SVM >= 0.90 ~ "Inv/Inv",
                            SVM <= 0.10 ~ "Std/Std",
                            SVM > 0.10 | SVM < 0.90 ~ "Std/Inv" ) ) %>%
  mutate(inv_freq = case_when(karyot == "Inv/Inv" ~ 0,
                              karyot == "Std/Std" ~ 2,
                              karyot == "Std/Inv" ~ 1) ) %>%
  mutate(geno = case_when(karyot == "Inv/Inv" ~ "Inv",
                          karyot == "Std/Std" ~ "Std",
                          karyot == "Std/Inv" ~ "Het") ) %>%
  mutate(Group = case_when(set == "CM" ~ "Alys",
                           set != "CM" ~ "Pris",) ) ->
  CM_predictions_parsed
####

CM_predictions_parsed %>%
  ggplot(
    aes(
      x=month(date_parsed),
      fill = geno
    )
  ) + 
  #geom_density(position = "fill") +
  geom_bar() +
  facet_wrap(~Group,  ncol = 1)->
  kar_bar

ggsave(kar_bar, file = "kar_bar.pdf")


####
CM_predictions_parsed %>% 
  group_by(date_parsed) %>%
  summarise(chr_tod = n()*2) -> tots


CM_predictions_parsed %>% 
  group_by(date_parsed,  Group) %>%
  summarise(gen_count = sum(inv_freq)) -> counts

left_join(tots, counts ) %>%
  mutate(prop = NA,
         prop_lo = NA,
         prop_h = NA) -> prop_table

for(i in 1:dim(prop_table)[1]){
  
  prop.test(prop_table$gen_count[i], 
            n = prop_table$chr_tod[i]) -> tmp
  
  prop_table$prop[i] = tmp$estimate
  prop_table$prop_lo[i] = tmp$conf.int[1]
  prop_table$prop_h[i] = tmp$conf.int[2]
  
}

prop_table %>%
  ggplot(aes(x=as.numeric(format(date_parsed, "%j")), y = prop, 
             #ymin =prop_lo, 
             ##ymax = prop_h
  )) +
  #geom_errorbar() +
  geom_line() +
  geom_point() +
  ggtitle("Frequency of the STD karyotypes in Alys/Pris Datasets") +
  xlab("Julian Day") +
  ylab("Frequency of STD") +
  facet_wrap(~Group,  ncol = 1)-> karyot_time

ggsave(karyot_time, file ="karyot_time.pdf", w = 5, h = 4)

#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data
#### ------> Investigate signal in the pool-seq data

### libraries
library(data.table)
library(SeqArray)
library(tidyverse)
library(car)

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

###########
###########
###########

#### load temperature data
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"
##-> weather.ave
##
### import the outlier SNPs
outSNPs <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")
#### load inversion markers
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
final_in2Lt_markers %<>%
  data.frame(markers = .) %>%
  separate(markers, into = c("chr", "pos", "type"), sep = "_")
final_in2Lt_markers$pos = as.numeric(final_in2Lt_markers$pos)

###
###

left_break = getData(chr="2L", start=2051609 , end=3096574)
right_break = getData(chr="2L", start=11584064 , end=13204668)

rbind(mutate(left_break, point = "left"), 
      mutate(right_break, point = "right")) %>%
  filter(locality == "VA_ch",
         pos %in% final_in2Lt_markers$pos) ->
  inv_bpoints

inv_bpoints %>%
  filter(year %in% 2016:2018) %>%
  ggplot(
    aes(
      x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=1-af_nEff,
      color=as.factor(variant.id)
    )) + 
  geom_line(alpha=0.5) +
  #geom_smooth(color = "black") +
  facet_grid(year~point, scales = "free_x") +
  scale_color_manual(values = rep("grey", length(unique(inv_bpoints$variant.id)))) +
  theme_classic() +
  theme(legend.position = "none") ->
  inv_breakpoints

ggsave(inv_breakpoints, file = "inv_breakpoints.pdf",  h=4, w=6)


### SIMULANS SIMULANS SIMULANS SIMULANS
### PLOT THE MEAN ALLELE FREQUENCY RELATIVE TO SIMULANS
padding=150000
starts= c(5155762, 6255762, 9505762)
ends= c(5255762, 6355762, 9605762)

win_5p = getData(chr="2L", start=5155762-padding , end=5255762+padding) %>%
  filter(#locality == "VA_ch",
         pos %in% outSNPs$pos)
win_6p = getData(chr="2L", start=6255762-padding , end=6355762+padding) %>%
  filter(#locality == "VA_ch",
         pos %in% outSNPs$pos)
win_9p = getData(chr="2L", start=9505762-padding , end=9605762+padding)%>%
  filter(#locality == "VA_ch",
         pos %in% outSNPs$pos)
### join datasets
rbind(mutate(win_5p, win = "1.win5"), 
      mutate(win_6p, win = "2.win6"),
      mutate(win_9p, win = "3.win9")) %>%   
  left_join(weather.ave) ->
  win_outliers_padding
### polarize with simulans
win_outliers_padding %>%
  filter( set %in% c("dgn"),
          sampleId == "SIM") %>%
  dplyr::select(chr, pos, sim_af = af) -> sim_af
#plotAllele frequencies
win_outliers_padding %>%
  left_join(sim_af) %>%
  mutate(af_neff_pol = case_when(
    sim_af == 0 ~ 1-af_nEff,
    sim_af == 1 ~ af_nEff)) %>%
  filter(year %in% 2016:2018,
         locality == "VA_ch"
         ) %>%
  group_by(pos,col,win) %>%
  summarize(mean_af_pol = mean(af_neff_pol, na.rm = T)) ->
  sum_Afs
###plot
ggplot() +
  geom_point(data=sum_Afs,
    aes(
      x=as.numeric(pos),
      y=1-mean_af_pol,
      color=col
    )
  ) +
  geom_vline(data=melt(data.frame(starts, ends, win = c("1.win5", "2.win6", "3.win9"))) ,
             aes(xintercept = value),
             color = "red",
             size = 0.8) +
  geom_point(size = 0.9) +
  ylab("Derived AF (rel. sim)") +
  theme_classic() +
  theme(legend.pos = "bottom") +
  facet_wrap(win~., scales = "free_x") -> af_sim_plot

ggsave(af_sim_plot, file = "af_sim_plot.pdf", w= 9, h = 3)


#### Trajectory vs. temperature
#### Trajectory vs. temperature
#### Trajectory vs. temperature
#### 
padding=0
starts= c(5155762, 6255762, 9505762)
ends= c(5255762, 6355762, 9605762)

win_5np = getData(chr="2L", start=5155762-padding , end=5255762+padding) %>%
  filter(#locality == "VA_ch",
    pos %in% outSNPs$pos)
win_6np = getData(chr="2L", start=6255762-padding , end=6355762+padding) %>%
  filter(#locality == "VA_ch",
    pos %in% outSNPs$pos)
win_9np = getData(chr="2L", start=9505762-padding , end=9605762+padding)%>%
  filter(#locality == "VA_ch",
    pos %in% outSNPs$pos)
### join datasets
rbind(mutate(win_5np, win = "1.win5"), 
      mutate(win_6np, win = "2.win6"),
      mutate(win_9np, win = "3.win9")) %>%   
  left_join(weather.ave) ->
  win_outliers_no_padding

##plot wins
win_outliers_no_padding %>%
  left_join(sim_af) %>%
  mutate(af_neff_pol = case_when(
    sim_af == 0 ~ 1-af_nEff,
    sim_af == 1 ~ af_nEff)) %>%
  filter(year %in% 2016:2018,
         locality == "VA_ch",
         col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant")
  ) %>% ggplot(
    aes(
      x=as.numeric(aveTemp/10),
      #x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=1-af_neff_pol,
      group=as.factor(variant.id),
      color=col,
    )) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(alpha=0.9, size = 0.5) +
  geom_line(alpha=0.5) +
  ylab("Derived AF (rel. sim)") +
  #geom_smooth(color = "black") +
  facet_grid(year~win, scales = "free_x") +
  #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
  #theme(legend.position = "none") +
  theme_bw() ->
  win_trajcs

ggsave(win_trajcs, file = "win_trajcs.pdf", h=4, w=8)
##plot wins





##plot SNPS
##plot SNPS
##plot SNPS

for(window in c("1.win5","2.win6","3.win9") ){
  
  win_outliers_no_padding %>%
    left_join(sim_af) %>%
    mutate(af_neff_pol = case_when(
      sim_af == 0 ~ 1-af_nEff,
      sim_af == 1 ~ af_nEff)) %>%
    filter(year %in% 2016:2018,
           win == window,
           #pos %in% filter(meand_derived_af, mean_pol_af > 0.5)$pos,
           locality == "VA_ch",
           col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant")
    ) %>% ggplot(
      aes(
        x=as.numeric(aveTemp/10),
        #x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
        y=(1-af_neff_pol),
        color=as.factor(year),
      )) + 
    geom_vline(xintercept = 20, linetype = "dashed") +
    geom_point(alpha=0.9, size = 0.5) +
    geom_smooth(se = F, 
                method = "lm"
                , color = "grey"
                ) +
    #geom_line(alpha=0.5) +
    #geom_smooth(color = "black") +
    facet_wrap(pos~gene, scales = "free") +
    ylab("Derived AF (rel. sim)") +
    #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
    #theme(legend.position = "none") +
    theme_bw() ->
    win_trajcs_snps
  
  ggsave(win_trajcs_snps, file = paste(window,"win_trajcs_snps.pdf", sep = "."), h=8, w=8)
  ##plot SNPS
}

#####boxplots
for(window in c("1.win5","2.win6","3.win9") ){
  
  win_outliers_no_padding %>%
    left_join(sim_af) %>%
    mutate(af_neff_pol = case_when(
      sim_af == 0 ~ 1-af_nEff,
      sim_af == 1 ~ af_nEff)) %>%
    filter(year %in% 2016:2018,
           win == window,
           #pos %in% filter(meand_derived_af, mean_pol_af > 0.5)$pos,
           locality == "VA_ch",
           col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant") ) %>% 
    mutate(temp_bin = RoundTo(aveTemp/10, 3 , "floor")) %>%
      ggplot(
      aes(
        x=as.factor(temp_bin),
         #x=as.numeric(aveTemp/10),
        #x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
        y=(1-af_neff_pol),
        #color=as.factor(year),
      )) + 
    geom_boxplot() +
    #geom_vline(xintercept = 20, linetype = "dashed") +
    #geom_point(alpha=0.9, size = 0.5) +
    #geom_smooth(se = F, 
    #            method = "lm"
    #            , color = "grey"
    #) +
    #geom_line(alpha=0.5) +
    #geom_smooth(color = "black") +
    facet_wrap(pos~gene, scales = "free") +
    ylab("Derived AF (rel. sim)") +
    #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
    #theme(legend.position = "none") +
    theme_bw() ->
    win_trajcs_snps.box
  
  ggsave(win_trajcs_snps.box, file = paste(window,"win_trajcs_snps.box.pdf", sep = "."), h=8, w=8)
  ##plot SNPS
}
#####boxplots


###### GEO
###### GEO
###### GEO
###### GEO
win_outliers_no_padding %>%
  left_join(sim_af) %>%
  mutate(af_neff_pol = case_when(
    sim_af == 0 ~ 1-af_nEff,
    sim_af == 1 ~ af_nEff)) %>%
  mutate(derived_maf = case_when(
    af_neff_pol > 0.5 ~ T,
    af_neff_pol <= 0.5 ~ F,
  )) %>%
  filter(
    lat > 26 & lat < 66,
    col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant"),
    set %in% c("DrosRTEC","DrosEU","CvilleSet")) %>%
  mutate(lat_bin = RoundTo(lat, 5, "floor") ) %>%
  ggplot(
    aes(
      x=as.factor(lat_bin),
      y=(1-af_neff_pol),
      fill=derived_maf
    )) +
    geom_boxplot() +
    ylab("Derived AF (rel. sim)") +
    facet_grid(.~win) ->
    global_mafs_wins
  ggsave(global_mafs_wins, file = "global_mafs_wins.pdf", h=5, w=8)
  
  ###### GEO TEMP
  ###### GEO TEMP
  ###### GEO TEMP
  ###### GEO TEMP
  ###### GEO TEMP
  ###### GEO TEMP
  
  win_outliers_no_padding %>%
    left_join(sim_af) %>%
    mutate(af_neff_pol = case_when(
      sim_af == 0 ~ 1-af_nEff,
      sim_af == 1 ~ af_nEff)) %>%
    mutate(derived_maf = case_when(
      af_neff_pol > 0.5 ~ T,
      af_neff_pol <= 0.5 ~ F,
    )) %>%
    filter(
      #lat > 26 & lat < 66,
      col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant"),
      set %in% c("DrosRTEC","DrosEU","CvilleSet")) %>%
    mutate(temp_bin = RoundTo(aveTemp/10, 5 , "floor")
           ) %>%
    .[complete.cases(.$temp_bin),] %>%
    ggplot(
      aes(
        x=as.factor(temp_bin),
        y=(1-af_neff_pol),
        fill=derived_maf
      )) +
    geom_boxplot() +
    ylab("Derived AF (rel. sim)") +
    facet_grid(.~win) ->
    global_mafs_wins.temp
  
  ggsave(global_mafs_wins.temp, file = "global_mafs_wins.temp.pdf", h=5, w=8)
  ###### GEO TEMP
  
  
  ###### GEO TEMP - comb
  ###### GEO TEMP - comb
  ###### GEO TEMP - comb
  ###### GEO TEMP - comb
  ###### GEO TEMP - comb
  ###### GEO TEMP - comb
  for(window in c("1.win5","2.win6","3.win9" ) ){
    
  win_outliers_no_padding %>%
    left_join(sim_af) %>%
    mutate(af_neff_pol = case_when(
      sim_af == 0 ~ 1-af_nEff,
      sim_af == 1 ~ af_nEff)) %>%
    mutate(derived_maf = case_when(
      af_neff_pol > 0.5 ~ T,
      af_neff_pol <= 0.5 ~ F,
    )) %>%
    filter(
      #lat > 26 & lat < 66,
      season %in% c("fall","spring"),
      win == window,
      col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant"),
      set %in% c("DrosRTEC","DrosEU"
                 #,"CvilleSet"
                 ),
      
      ) %>%
    mutate(temp_bin = RoundTo(aveTemp/10, 5 , "floor"),
           lat_bin = RoundTo(lat, 5, "floor")) %>%
    #.[complete.cases(.$temp_bin),] %>%
    ggplot(
      aes(
        x=as.factor(lat_bin),
        y=1-af_neff_pol,
        fill=as.factor(season)
      )) +
    geom_boxplot() +
      ggtitle(window, subtitle = "DEST only")+
    ylab("Derived AF (rel. sim)") +
      facet_wrap(pos~gene, scales = "free") ->
      global_mafs_wins.temp.lat
  
  ggsave(global_mafs_wins.temp.lat, file = paste(window, "global_mafs_wins.temp.lat.pdf", sep = "."), h=8, w=8)
  }
 
  
  
  
  
  
  
  
  #OLD 
  #  #OLD   #OLD 
  #    #OLD 
  #      #OLD 
######

#ex = getData(chr="2L", start=9538630 , end=9538630)

  for(window in c("1.win5","2.win6","3.win9" ) ){
    
  win_outliers_no_padding %>%
  left_join(sim_af) %>%
  mutate(af_neff_pol = case_when(
    sim_af == 0 ~ 1-af_nEff,
    sim_af == 1 ~ af_nEff)) %>%
  filter(
         lat > 26 & lat < 66,
         win == window,
         #long < 20 & long > -90,
         col %in% c("missense_variant","3_prime_UTR_variant", "5_prime_UTR_variant"),
         set %in% c("DrosRTEC","DrosEU"),
         season %in% c( "fall", "spring")) %>%
  #mutate( lat_binned = RoundTo(lat, 5, "floor")) %>%
  ggplot(
    aes(
      x=as.numeric(lat),
      #x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=(1-af_neff_pol),
      color=as.factor(continent),
      shape = season,
      linetype = season
      
    )) +
      geom_point(alpha=2, size = 0.5) +
      geom_smooth(se = F, method = "lm") +
      #geom_line(alpha=0.5) +
      facet_wrap(pos~gene, scales = "free") +
      #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
      #theme(legend.position = "none") +
      theme_bw() ->  ex_geo

ggsave(ex_geo, file = paste(window, "ex_geo.pdf", sep = "."), h=5, w=8)
  }
  

#### LD investigations
#### LD investigations
#### LD investigations
#### LD investigations
#### LD investigations
#### LD investigations
#### 
#### 
### load libraries
library(patchwork)
library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(data.table)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)

load("../7.LD/merged.ld.Rdata")
output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")
outSNPs <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")


chr="2L"
for(pos in c(6632309, 10184864)){
  message(pos)
  
  ld_df %>%
    filter(CHR_A == chr,
           CHR_B == chr,
           BP_A == pos | BP_B == pos,
           R2 > 0.6
    ) %>%
    mutate(Type = "all" ) %>%
    .[order(.$BP_A)] %>%
    mutate(pair_id = 1:dim(.)[1],
           bp_dist = abs(BP_A-BP_B)) %>%
    .[order(.$BP_A, .$bp_dist)] -> gen_ld_all
  
  gen_ld_all %>%
    left_join(., dplyr::select(outSNPs, BP_A=pos, AnnotA=Annotation, geneA=Gene_Name)) %>%
    left_join(., dplyr::select(outSNPs, BP_B=pos, AnnotB=Annotation, geneB=Gene_Name)) %>%
    filter(AnnotA == "missense_variant", AnnotB == "missense_variant" ) %>%
    mutate(Type = "NS_vars" ) %>%
    .[order(.$BP_A)] %>%
    mutate(pair_id = 1:dim(.)[1],
           bp_dist = abs(BP_A-BP_B)) %>%
    .[order(.$BP_A, .$bp_dist)] -> gen_ld_NS
  

  rbind(gen_ld_all, gen_ld_NS, fill = T) %>%
    ggplot(
      aes(
        x=BP_A,
        xend=BP_B,
        y=pair_id,
        yend=pair_id,
        color=R2
      )
    ) +
    geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
    geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
    scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.80) +
    geom_segment() +
    facet_wrap(~Type, scale = "free_y") ->
    gen_plot
  
  ggsave(gen_plot, file = paste(pos, "gen_plot.pdf", sep = "") )
  
  print(unique(c(gen_ld_NS$geneA, gen_ld_NS$geneB)))

}

  

