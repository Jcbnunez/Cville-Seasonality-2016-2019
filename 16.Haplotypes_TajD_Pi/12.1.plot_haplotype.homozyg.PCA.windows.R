library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)
library(foreach)
library(SeqArray)
library(doMC)
library(car)
library(DescTools)
library(ape)
library(ggtree)
library(aplot)
library(forcats)
registerDoMC(2)

###
# setwd("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/")
# load the inversion markers
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")

#load("/scratch/yey2sn/Overwintering_ms/4.GML_plots/PEAKS_for_ANALYSIS.Rdata")
#PEAKS_for_ANALYSIS

#outlier_haplowins = 
#  data.frame(win.name = c("win_4.6", "win_5.1", "win_6.2", "win_6.8", "win_9.5" ),
#             start = c(4656622, 5105919, 6155931, 6805798, 9505855 ),
#             end = c(4805715, 5255685, 6355509, 6905746, 9605419))

final.windows.pos = 
  data.frame(win.name = c("left", "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6", "right" ),
             mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )



# load outlier GLMs
#glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"
##load model files
#mods_fin <- fread("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/final_models.txt")
#mods_fin$mod_var -> models

#load("/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/real.dat.models.Rdata")
#real.dat.models

# load best Cville model
file.mod <- "/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata"
glm.out <- get(load(file.mod))

glm.out %>%
  filter(perm == 0) %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  filter(chr == "2L",
         rnp < 0.05) %>% 
  mutate(win = case_when(
    pos > 2800000 & pos < 3200000 ~ "win_3.1",
    pos > 4470000 & pos < 4870000 ~ "win_4.7",
    pos > 4920000 & pos < 5320000 ~ "win_5.1",
    pos > 6000000 & pos < 6400000 ~ "win_6.1",
    pos > 6600000 & pos < 7000000 ~ "win_6.8",
    pos > 9400000 & pos < 9800000 ~ "win_9.6",
    pos < 3e6 & SNP_id %in% final_in2Lt_markers ~ "left",
    pos > 11e6 & SNP_id %in% final_in2Lt_markers ~ "right")) -> 
  glm.outliers.2L.wins

glm.outliers.2L.wins %>%
  filter(win %in% c("win_3.1", "win_4.7" ,"win_5.1", "win_6.1" , "win_6.8" ,"win_9.6", "left", "right"))  -> 
  glm.outliers.2L.wins.flt

#outliers_glm %<>% separate(V1, into = c("chr", "pos", "type"))

#load windows data
#load("./AF_dat_summ_id.dat.Rdata")
#AF_dat_summ_id$win %>% table
#SIM_AF %>% head
###

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM.toPredict_INV_STD.new.Apr1.2022.Rdata")

sample_metadata %>%
  filter(pop != "NY") %>%
  filter(INV_STATUS %in% c("INV", "STD")) %>%
  dplyr::select(hap_name=sampleId, karyo=INV_STATUS, pop ) ->
  Nat_pops_ids
  
STD_DGRP = fread("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/STD_DGRP_OnlyNames.txt", header = F)
STD_DGRP %<>% 
  dplyr::select(hap_name=V1) %>%
  mutate(karyo = "STD", pop = "DGRP")

INV_DGRP = fread("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/INV_DGRP_OnlyNames.txt", header = F)
INV_DGRP %<>% 
  dplyr::select(hap_name=V1) %>%
  mutate(karyo = "INV", pop = "DGRP")

rbind(Nat_pops_ids, STD_DGRP, INV_DGRP) -> samps_info_pop


####
####
vcf_file <- "./MC.DGRP.Taylor.merged.readyForImport.recode.vcf.gz"
vcf <- read.vcfR( vcf_file, verbose = TRUE )

my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE)
##my_genind <- vcfR2genind(vcf)
#---> AF_dat_summ_id

##left = filter(AF_dat_summ_id, win == "left", SNP_id %in% final_in2Lt_markers)
##right = filter(AF_dat_summ_id, win == "right", SNP_id %in% final_in2Lt_markers)
##win5 = filter(AF_dat_summ_id, win == "win_5", SNP_id %in% glm.outliers.2L$SNP_id)
##win6 = filter(AF_dat_summ_id, win == "win_6", SNP_id %in% glm.outliers.2L$SNP_id)
##win9 = filter(AF_dat_summ_id, win == "win_9", SNP_id %in% glm.outliers.2L$SNP_id)



### left
left_snps_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
  filter(glm.outliers.2L.wins.flt, win == "left")$SNP_id)]) 

markers_left = DNAbin2genind(x = 
                            my_dnabin1[,which(colnames(my_dnabin1) %in% 
                            filter(glm.outliers.2L.wins.flt, win == "left")$SNP_id)]  
                            ,polyThres = 0.00)
locNames(markers_left) = left_snps_names

### right
right_snps_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
  filter(glm.outliers.2L.wins.flt, win == "right")$SNP_id)]) 

markers_right = DNAbin2genind(x = 
my_dnabin1[,which(colnames(my_dnabin1) %in% 
filter(glm.outliers.2L.wins.flt, win == "right")$SNP_id)]  
,polyThres = 0.00)

locNames(markers_right) = right_snps_names


## win_3.1
win_3.1_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_3.1")$SNP_id)]) 

markers_win_3.1 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                                        filter(glm.outliers.2L.wins.flt, 
                                                               win == "win_3.1")$SNP_id)]  
                                ,polyThres = 0.00)

locNames(markers_win_3.1) = win_3.1_names

###

## win_4.7
win_4.7_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_4.7")$SNP_id)]) 

markers_win_4.7 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                                        filter(glm.outliers.2L.wins.flt, 
                                                               win == "win_4.7")$SNP_id)]  
                                ,polyThres = 0.00)

locNames(markers_win_4.7) = win_4.7_names

###
## win_5.1
win_5.1_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_5.1")$SNP_id)]) 

markers_win_5.1 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                                    filter(glm.outliers.2L.wins.flt, 
                                                           win == "win_5.1")$SNP_id)]  
                                ,polyThres = 0.00)

locNames(markers_win_5.1) = win_5.1_names

###

## win_6.1
win_6.1_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_6.1")$SNP_id)]) 

markers_win_6.1 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                                        filter(glm.outliers.2L.wins.flt, 
                                                               win == "win_6.1")$SNP_id)]  
                                ,polyThres = 0.00)

locNames(markers_win_6.1) = win_6.1_names

## win_6.8
win_6.8_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_6.8")$SNP_id)]) 

markers_win_6.8 = DNAbin2genind(x = my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                                        filter(glm.outliers.2L.wins.flt, 
                                                               win == "win_6.8")$SNP_id)]  
                                ,polyThres = 0.00)

locNames(markers_win_6.8) = win_6.8_names


## win_9.6
win_9.6_names =colnames(
  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                      filter(glm.outliers.2L.wins.flt, win == "win_9.6")$SNP_id)]) 

markers_win_9.6 = DNAbin2genind(x = 
                                  my_dnabin1[,which(colnames(my_dnabin1) %in% 
                                  filter(glm.outliers.2L.wins.flt, win == "win_9.6")$SNP_id)]  
                                  ,polyThres = 0.00)

locNames(markers_win_9.6) = win_9.6_names



### Join datasets
datasets = list(left_w = markers_left,
                right_w = markers_right,
                win_3.1 = markers_win_3.1,
                win_4.7 = markers_win_4.7,
                win_5.1 = markers_win_5.1,
                win_6.1 = markers_win_6.1,
                win_6.8 = markers_win_6.8,
                win_9.6 = markers_win_9.6
                )

#################################################
##### PLOT JOINNT FIGURE
joint_figure = foreach(i = 1:length(datasets), .combine = "rbind" )%dopar%{
  #for(i in 1:length(datasets)){
  data_in = datasets[[i]]
  
  message(names(datasets)[i])
  
  data_in@tab %>%
    as.data.frame %>% 
    as.data.frame -> data_in_tab
  
  data_in_tab[,seq(from=1, to=dim(data_in_tab)[2], by=2)] -> data_in_tab_loc
  dim(data_in_tab_loc)
  
  colnames(data_in_tab_loc) %>%
    data.frame(samp_hap = .) %>% 
    mutate(loci_id = 1:dim(.)[1]) ->  metadat_loci
  
  data_in_tab_loc %>%
    mutate(samp_id = rownames(.)) %>%
    melt(id = c("samp_id"), variable.name = "loci_id") ->
    data_in_tab_loc_melt
  
  data_in_tab_loc_melt$loci_id = as.integer(data_in_tab_loc_melt$loci_id )
  
  data_in_tab_loc_melt %<>%
    mutate(hap_name = gsub("_0$|_1$", "", data_in_tab_loc_melt$samp_id))
  
  ##add metadat
  data_in_tab_loc_melt %<>%
    left_join(samps_info_pop) %>%
    .[complete.cases(.),]
  
  left_join(data_in_tab_loc_melt, metadat_loci) %>% 
    separate(samp_hap, remove = F, into = c("SNP_id", "base"), sep = "\\.") %>%
    mutate(win =  names(datasets)[i]) -> objt
  
  return(objt)
  
}

##prepare plot
## polarization step  

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

### run AF collector

rbind(
  mutate(getData(chr="2L", start=2051609 , end=3096574), win = "left") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           SNP_id %in% final_in2Lt_markers),
  
  mutate(getData(chr="2L", start=11584064 , end=13204668), win = "right") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           SNP_id %in% final_in2Lt_markers),
  
  getData(chr="2L", start=2800000, end=3200000) %>% mutate(win="win_3.1") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos),
  
  getData(chr="2L", start=4470000, end=4870000) %>% mutate(win="win_4.7") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos),
  
  getData(chr="2L", start=4920000, end=5320000) %>% mutate(win="win_5.1")%>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos),
  
  getData(chr="2L", start=6000000, end=6400000) %>% mutate(win="win_6.1") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos),

  getData(chr="2L", start=6600000, end=7000000) %>% mutate(win="win_6.8") %>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos),
    
  getData(chr="2L", start=9400000, end=9800000) %>% mutate(win="win_9.6")%>%
    mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
    filter(sampleId == "SIM",
           pos %in% glm.outliers.2L.wins.flt$pos)
  ) %>%
  arrange(pos) %>%
  mutate(tidy_annot = case_when(
    col %in% c("3_prime_UTR_variant","5_prime_UTR_variant") ~ "UTR",
    col %in% c("downstream_gene_variant","intergenic_region", "upstream_gene_variant") ~ "intergenic",
    col %in% c("intron_variant","splice_region_variant&intron_variant") ~ "intron",
    col %in% c("missense_variant") ~ "NS",
    col %in% c("synonymous_variant", "splice_region_variant&synonymous_variant") ~ "S",
  )) ->
  sim_polarity

sim_polarity %>% 
  dplyr::select(sampleId, variant.id, gene, SNP_id, tidy_annot, af) ->
sim_polarity.flt

### Merge datasets
joint_figure %>% 
  left_join(sim_polarity.flt, by = c("SNP_id") ) %>% 
  mutate(value_pol = case_when(
    af == 0 ~ 1-as.numeric(value),
    af == 1 ~ as.numeric(value))) %>% 
  mutate(polarity = case_when(
    value_pol == 0 ~ "Ancestral",
    value_pol == 1 ~ "Derived")) ->
  joint_figure_polarized

joint_figure_polarized$win = factor(joint_figure_polarized$win, 
                                    levels = c("left_w",
                                               "win_3.1",
                                               "win_4.7",
                                               "win_5.1",
                                               "win_6.1",
                                               "win_6.8",
                                               "win_9.6",
                                               "right_w")
                                    )


save(joint_figure_polarized, file = "joint_figure_polarized.Rdata")
load("./joint_figure_polarized.Rdata")

### filter by num of sites
select_pops = c("CM","PA","ME","DGRP","France","Netherlands")

#Cameroon       Carib          CM        DGRP      France      Guinea 
#20232        8816      388596      499818       22780       12642 
#ME Netherlands          PA      Zambia 
#190162       37608      265974       72936 

joint_figure_polarized %>%
  filter(pop %in% select_pops) %>% 
  group_by(SNP_id) %>% 
  summarize(ALL= n()) -> allSNPS

joint_figure_polarized %>%
  filter(is.na(value_pol), pop %in% select_pops) %>% 
  group_by(SNP_id) %>% 
  summarize(Miss= n()) -> missSNPS

full_join(allSNPS, missSNPS) -> miss_analysis

miss_analysis$Miss[is.na(miss_analysis$Miss)] = 0 

miss_analysis %<>%
  mutate(perc_miss = as.numeric(Miss)/as.numeric(ALL)) 

 miss_analysis$ALL %>% quantile(0.7) -> tresh_all
 miss_analysis$perc_miss %>% quantile(0.7) -> thres_miss
 
 miss_analysis %>%
  filter(ALL >= tresh_all,
         perc_miss <= thres_miss) %>% 
  .$SNP_id -> id_id_wins

### filter by inds
### 

joint_figure_polarized %>%
  filter(pop %in% select_pops) %>% 
  group_by(hap_name) %>% 
  summarize(ALL= n()) -> inds_snps_samps_density

quantile(inds_snps_samps_density$ALL)

inds_snps_samps_density %>%
  filter(ALL >= 5440) %>% 
  .$hap_name -> keep_samps


#######
####### Make filterd object
joint_figure_polarized %>%
filter(SNP_id %in% id_id_wins,
       hap_name %in% keep_samps ) ->
  joint_figure_polarized_missingdat_clean


#######
#######
#######
#######
#######
#######

joint_figure_polarized_missingdat_clean %>%
  filter(pop %in% select_pops) %>% 
  .[complete.cases(.$value_pol),] %>%
ggplot(
    aes(
      x=as.factor(loci_id),
      y=samp_id,
      fill = polarity
    )
  ) + geom_tile(size = 0.1) +
  facet_grid(karyo+pop~win, scales = "free",  
             space = "free" 
             #ncol = 1, shrink = F
             ) +
  ggtitle("joint hapblocks") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  ->
  karyo_plot_joint

ggsave(karyo_plot_joint, file =  "karyo_plot_joint.2.pdf", w = 12, h = 8)

############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering
############# PHYLOGENETICS and clustering

samp_names =
  data.frame(rbind(
cbind(samp=paste(filter(samps_info_pop, karyo == "STD", 
                        hap_name %in% keep_samps,
                        pop %in% select_pops)$hap_name, 
                 "0", sep = "_"), type = "STD"),
cbind(samp=paste(filter(samps_info_pop, karyo == "INV", 
                        hap_name %in% keep_samps,
                        pop %in% select_pops)$hap_name, 
                 "0", sep = "_"), type = "INV"),
cbind(samp=paste(filter(samps_info_pop, karyo == "STD", 
                        hap_name %in% keep_samps,
                        pop %in% select_pops)$hap_name, 
                 "1", sep = "_"), type = "STD"),
cbind(samp=paste(filter(samps_info_pop, karyo == "INV", 
                        hap_name %in% keep_samps,
                        pop %in% select_pops)$hap_name, 
                 "1", sep = "_"), type = "INV")
  ))


set.seed(123456)

select_haps = data.frame( 
  rbind(filter(samp_names, type == "INV"), 
        sample_n(filter(samp_names, type == "STD" ), 30, replace = FALSE))
)

joint_figure_polarized_missingdat_clean %>%
  filter(samp_id %in% select_haps$samp) %>% 
  dcast(samp_id~SNP_id, value.var = "value_pol")  ->
  windws_snp_matrix_clean

### Parts of Figure 4
save(
  sim_polarity,
  samp_names,
  windws_snp_matrix_clean,
  joint_figure_polarized_missingdat_clean,
  #file = "/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata"
  file = "./haplowins_pt1.Rdata"
  
)

#load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/haplowins_pt1.Rdata")
load("./haplowins_pt1.Rdata")
## add rownames
rownames(windws_snp_matrix_clean) = windws_snp_matrix_clean$samp_id

analyses_types = list(
  all=c("left",
        "win_3.1",
        "win_4.7",
        "win_5.1",
        "win_6.1",
        "win_6.6",
        "win_9.6",
        "right") #,
  #win3.1=c("win_3.1"),
  #win5.1=c("win_5.1"),
  #win9.6=c("win_9.6")
  )

foreach(i=1:length(analyses_types))%do%{

  message(names(analyses_types)[i])
  
analysis = names(analyses_types)[i]
target_snps <- filter(sim_polarity, win %in% as.character(analyses_types[[i]]) )$SNP_id 

data_in <- windws_snp_matrix_clean[, which(colnames(windws_snp_matrix_clean) %in% target_snps)  ] 

actual_sel_snps <- colnames(data_in)

# dplyr::select(windws_snp_matrix_clean, !(samp_id), filter(SNP_guide_metadata, window == "win5" )$SNP_id )
D_all <- dist( data_in )
tre_all <- njs(D_all)

tree_all_plot <- ggtree(tre_all, ignore.negative.edge=TRUE) + 
  #geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

tree_all_plot <- tree_all_plot %<+% joint_figure_polarized_missingdat_clean + geom_tippoint(aes(color=pop))


is_tip <- tre_all$edge[,2] <= length(tre_all$tip.label)
ordered_tips <- tre_all$edge[is_tip, 2]
tre_all$tip.label[ordered_tips]  -> tree_order

joint_figure_polarized_missingdat_clean %>%
  #filter(samp_id %in% select_haps$samp[grep("line", select_haps$samp,  invert = T)] ) %>% 
  mutate(samp_id_fct = factor(samp_id, levels = tree_order)) %>% 
  .[complete.cases(.$samp_id_fct),] %>%
  ggplot(
    aes(
      x=as.factor(loci_id),
      y=samp_id_fct,
      fill = polarity
    )
  ) + geom_tile(size = 0.1) +
  facet_grid(.~win, scales = "free", space = "free"
             #ncol = 1, shrink = F
  ) +
  ggtitle(analysis) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  ->
  karyo_plot_joint_all

karyo_plot_joint_all %>% 
  insert_left(tree_all_plot, width = 0.1)  -> hap_tree_plots_all

ggsave(hap_tree_plots_all, file = paste(analysis, "hap_tree_plots.all.pdf", sep = "."), h = 6, w = 10)

}

#### Add annotation 
#### 

target_snps <- filter(sim_polarity, win %in% analyses_types[[1]] )$SNP_id 
data_in <- windws_snp_matrix_clean[, which(colnames(windws_snp_matrix_clean) %in% target_snps)  ]
actual_sel_snps <- colnames(data_in)

joint_figure_polarized_missingdat_clean %>%
  filter(SNP_id %in% actual_sel_snps) %>%
  group_by(SNP_id) %>%
  slice_head() %>%
  ggplot(aes(
    x=as.factor(loci_id),
    y=1,
    fill = tidy_annot
  )) + geom_tile(size = 0.1) +
  facet_grid(.~win, scales = "free", space = "free"
  ) +
  #ggtitle(analysis) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")  ->
  annots_plot
ggsave(annots_plot, file =  "annots_plot.pdf",  h = 1.5, w = 10)

### annotate haplotags
### ### annotate haplotags
### annotate haplotags
### annotate haplotags
### annotate haplotags
### annotate haplotags
### annotate haplotags

target_snps <- filter(sim_polarity, win %in% analyses_types[[1]] )$SNP_id 
data_in <- windws_snp_matrix_clean[, which(colnames(windws_snp_matrix_clean) %in% target_snps)  ]
actual_sel_snps <- colnames(data_in)

### This needs to complete the haplotags code
### ### This needs to complete the haplotags code
### This needs to complete the haplotags code
### This needs to complete the haplotags code
### This needs to complete the haplotags code
### This needs to complete the haplotags code
### This needs to complete the haplotags code
### This needs to complete the haplotags code

#load("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/haplo_tags_SNPids.Rdata")
load("./haplo_tags_SNPids.Rdata")


haplo_tags_SNPids

joint_figure_polarized_missingdat_clean %>%
  mutate(haplotag = case_when(
    joint_figure_polarized_missingdat_clean$SNP_id %in% haplo_tags_SNPids$SNP_id ~ "yes",
    !(joint_figure_polarized_missingdat_clean$SNP_id %in% haplo_tags_SNPids$SNP_id) ~ "no"
  )) -> hap_tags

hap_tags %>%
  filter(SNP_id %in% actual_sel_snps) %>%
  group_by(SNP_id) %>%
  slice_head() %>%
  ggplot(aes(
    x=as.factor(loci_id),
    y=1,
    fill = haplotag
  )) + geom_tile(size = 0.1) +
  facet_grid(.~win, scales = "free", space = "free"
  ) +
  #ggtitle(analysis) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")  ->
  haplotag_plot
ggsave(haplotag_plot, file =  "haplotag_plot.pdf",  h = 1.5, w = 10)

