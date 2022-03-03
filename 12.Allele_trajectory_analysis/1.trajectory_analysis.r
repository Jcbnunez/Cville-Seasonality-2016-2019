### Genetics and Environemtns characterization of Genetic targets

## Part 1 -- Identify inversion status

library(tidyverse)
library(magrittr)
library(data.table)

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

#min_pos > 4.7e6 & max_pos < 5.3e6 ~ "1.win5",
#min_pos > 5.8e6 & max_pos < 7e6 ~ "2.win6",
#min_pos > 9.0e6 & max_pos < 10.6e6 ~ "3.win10")) 

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



###
###
win_5 = getData(chr="2L", start=4.7e6 , end=5.3e6) %>%
  filter(locality == "VA_ch",
         pos %in% outSNPs$pos)
win_6 = getData(chr="2L", start=5.8e6 , end=7e6) %>%
  filter(locality == "VA_ch",
         pos %in% outSNPs$pos)
win_10 = getData(chr="2L", start=9.0e6 , end=10.6e6)%>%
  filter(locality == "VA_ch",
         pos %in% outSNPs$pos)


###
###
rbind(mutate(win_5, win = "1.win_5"), 
      mutate(win_6, win = "2.win_6"),
      mutate(win_10, win = "3.win_10")) %>%   
  left_join(weather.ave) ->
  win_outliers

win_outliers %>%
  filter(year %in% 2016:2018,
         col %in% c("synonymous_variant",
                    "intron_variant",
                    "missense_variant")
         ) %>%
  ggplot(
    aes(
      #x=as.numeric(aveTemp/10),
      x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=1-af_nEff,
      group=as.factor(variant.id),
      color=col
    )) + 
  geom_line(alpha=0.5) +
  #geom_smooth(color = "black") +
  facet_grid(year~win, scales = "free_x") +
  #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
  #theme(legend.position = "none") +
  theme_classic() ->
  win_trajcs

ggsave(win_trajcs, file = "win_trajcs.pdf", h=4, w=8)

win_outliers%>%
  filter(gene %in% "Sur",
         year %in% 2016:2018,
         col == "missense_variant",
         af_nEff > 0.7) 
  

win_outliers%>%
  filter(gene %in% "Sur",
         year %in% 2016:2018) %>%
  ggplot(
    aes(
      #x=as.numeric(aveTemp/10),
      x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=1-af_nEff,
      group=as.factor(variant.id),
      color=col
    )) + 
  geom_line(alpha=0.5) +
  #geom_point()+
  #geom_smooth(color = "black") +
  facet_grid(year~pos, scales = "free_x") +
  #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
  theme_classic()  +
  theme(legend.position = "top",
        text = element_text(size = 8)) ->
  win_trajcs_high

ggsave(win_trajcs_high, file = "win_trajcs_high.pdf", h=4, w=9)


###### Taj D zooms

TajD6 = getData(chr="2L", start=6.07e6 , end=6.69e6) %>%
  filter(locality == "VA_ch",
         pos %in% outSNPs$pos)
TajD10 = getData(chr="2L", start=9.64e6 , end=10.0e6) %>%
  filter(locality == "VA_ch",
         pos %in% outSNPs$pos)

rbind(mutate(TajD6, win = "1.TajD6"), 
      mutate(TajD10, win = "2.TajD10")) %>%   
  left_join(weather.ave) ->
  Taj_outliers


Taj_outliers %>%
  .[grep("^CG" ,.$gene, invert = T)] %>%
  filter(year %in% 2016:2018,
         col %in% c(#"synonymous_variant",
                   # "intron_variant",
                    "missense_variant"
                    )) %>%
  ggplot(
    aes(
      #x=as.numeric(aveTemp/10),
      x=as.numeric(format(as.Date(collectionDate, format = c("%m/%d/%Y")), "%j")),
      y=1-af_nEff,
      color=as.factor(variant.id),
      linetype=as.factor(year)
    )) + 
  geom_line(alpha=0.5) +
  #geom_smooth(color = "black") +
  facet_wrap(~gene, scales = "free_x") +
  #scale_color_manual(values = rep("grey", length(unique(win_outliers$variant.id)))) +
  theme_classic() +
  theme(legend.position = "none")  ->
  Taj_outliers_win_trajcs

ggsave(Taj_outliers_win_trajcs, file = "Taj_outliers_win_trajcs.pdf", h=4, w=8)

#####

ex = getData(chr="2L", start=6469925 , end=6469925)

ex %>%
  filter(long < -50 & long > -93,
         set %in% c("DrosRTEC","DrosEU","CvilleSet")) %>%
  ggplot(aes(
    x=lat,
    y=af_nEff,
    color=season
  )) +
  geom_smooth(method = "lm", color = "grey") +
  geom_point() ->
  ex_geo

ggsave(ex_geo, file = "ex_geo.pdf")

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

  

