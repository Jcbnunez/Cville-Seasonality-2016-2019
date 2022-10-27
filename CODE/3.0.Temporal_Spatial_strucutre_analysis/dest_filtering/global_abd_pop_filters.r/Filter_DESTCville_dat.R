#Code to generate worldwide and population specific filters for the DEST+Cville data
#

library(data.table)
library(SeqArray)
library(tidyverse)
library(lubridate)
library(magrittr)
library(gmodels)

#Part 1. Download dataset
# load recombination rates metadata
areas_of_low_rec <- fread("/project/berglandlab/Dmel_genomic_resources/Recombination_Rates/hglift_Dmelv6_areas_with0recomb.bed")
names(areas_of_low_rec) = c("chr","begin","end","Rate_cM_Mb")


# load  sample metadata
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

#Make corrections to the metadata --- some names are not standardized
samps$city[grep("Odesa", samps$city )] = "Odessa"
samps$city[grep("Charlotttesville", samps$city )] = "Charlottesville"
samps$city[grep("Yesiloz", samps$city )] = "Yesiloz"
samps$city[grep("Chornobyl", samps$city )] = "Chernobyl"
samps$city[grep("Kyiv", samps$city )] = "Kyiv"
samps$city[grep("Uman", samps$city )] = "Uman"
samps$city[grep("valday", samps$city )] = "Valday"

#Select sets
#This will exclude all DGN samples and other insilico pools
samps <- rbind(samps[set=="DrosRTEC"],
               samps[set=="DrosEU"],
               samps[set=="CvilleSet"]
)

#Make some modification to the collection metadata
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps[,Date:=date(paste(year, month, day, sep="-"))]

#Import the GDS object into R
### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", 
                    allow.duplicate=T)

### obtain the metadata of all common SNPs... to be saved into object SNP.dt
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id")) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

#characterize multialleles
snp.dt %>%
  group_by(nAlleles) %>%
  summarize(N= n())

#Filtering only biallelic SNPs
snp.dt <- snp.dt[nAlleles==2]
#This step filters the geno file based on the SNPid from snp.dt
seqSetFilter(genofile, snp.dt$id)

#How many SNPs are found within chromosomes
snp.dt %>%
  group_by(chr) %>%
  summarize(N= n())

#add extra metadata to the snp.dt object
#Add the mean AF for all SNPs. A global AF filter
snp.dt[,mean.af:=seqGetData(genofile, "annotation/info/AF")$data]
#snp.dt <- snp.dt[mean.af>.1 & mean.af<.9]

#Describe allelelic states
snp.dt$mean.af %>% quantile %>% round(3)

snp.dt %>%
  ggplot(
    aes(mean.af)
  ) +
  geom_histogram() +
  facet_wrap(~chr,
             scales = "free",
             nrow = 1) +
  xlab("Mean AF") ->
  chrs_AF

ggsave(chrs_AF,
       file = "chrs_AF.png",
       width = 12,
       height = 2.5)

# make the chromosome column a guiding column
setkey(snp.dt, chr)
#Extract all autosomes but exclude "X".
snp.dt <- snp.dt[J(c("2L", "2R", "3L", "3R"))]

#####################################
#Evaluate areas of low recombination ---> can jump to saved object!
#####################################
chrs = snp.dt$chr %>% unique

#make list
snps_in_low_rec_list = list()

#loop over chromosomes
for(i in 1:length(chrs)){
  
range <-  data.table(start = areas_of_low_rec$begin[which(areas_of_low_rec$chr == chrs[i])],
                     end = areas_of_low_rec$end[which(areas_of_low_rec$chr == chrs[i])]) 

#subset data ... for spped
dat <- snp.dt[which(snp.dt$chr == chrs[i])]

#used a vectorized solution for speed
inrange=dat$SNP_id[dat$pos %inrange% range]

snps_in_low_rec_list[[i]] = inrange

}

#unlist
snps_in_low_rec_vector = unlist(snps_in_low_rec_list)

#add data to the snp.dt object
snp.dt %<>%
  mutate(Low_rec_zone = ifelse(.$SNP_id %in% 
                                 snps_in_low_rec_vector,
                                 TRUE, FALSE)) 

snp.dt %>%
  group_by(chr, Low_rec_zone) %>%
  summarise(Mean = quantile(2*mean.af*(1-mean.af), 0.5),
            q75= quantile(2*mean.af*(1-mean.af), 0.75),
            q25= quantile(2*mean.af*(1-mean.af), 0.25)
            #Low = ci(mean.af, confidence=0.9)[2],
            #High = ci(mean.af, confidence=0.9)[3],
            ) %>% 
  ggplot(aes(
    x=chr,
    y=Mean,
    ymin=q25,
    ymax=q75,
    fill=Low_rec_zone
  )) + 
  geom_errorbar(aes(color = Low_rec_zone),
                width = 0.4,
                position = position_dodge(width = 0.5))  +
  scale_fill_manual(name = "Low rec. region", values = c("blue","red")) +
  scale_color_manual(name = "Low rec. region", values = c("blue","red")) +
  geom_point(shape = 21,
             size = 4,
             position = position_dodge(width = 0.5)) +
  ylab(expression( H[italic(E)] )) +
  xlab("Chromosome Arm")->
  rec_af

ggsave(rec_af,
       file = "rec_af.pdf",
       width = 5,
       height = 4)

# Save Intermediary work:
#save(snp.dt, file = "snp.dt.savePoint.Rdata")  

### snp.dt <--------- contains SNP metadata after rep abd recon filtering

## Part 2. Filter populations based on EC data

### select sites
seqSetFilter(genofile, 
             sample.id=samps$sampleId,
             snp.dt$id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad/dp #if poolSNP
dim(dat)  

## Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

## Calculate effective coverage
### Generate a function to estimate effective coverage
### Nc = (1/N + 1/R) - 1
calc_effcov = function(NomCOV, FlyPool) {
  
  EC1 = ((2*FlyPool)*NomCOV)
  EC2 = ((2*FlyPool)+NomCOV)
  Nc=(EC1/EC2)
  return(Nc)
  
}

#######
dp %>%
  .[,sample(dim(.)[2], 100000 )] %>%
  as.data.table() %>%
  mutate(sampleId = rownames(dp)) %>%
  left_join(samps[,c("sampleId","nFlies")]) %>%
  melt(id = c("sampleId","nFlies")) %>% 
  mutate(EFFCOV = calc_effcov(value, nFlies) ) %>% 
  group_by(sampleId) %>%
  summarise(MeanEC = mean(EFFCOV, na.rm = T)) -> EFFCOV_samps

#### merge Neff with samps
samps %<>%
  left_join(EFFCOV_samps)

#############################
#### add sample filtering
samps %<>%
  mutate(Cville_only_filter  = ifelse(.$set == "CvilleSet" & .$MeanEC < 30, "fail", "pass"),
         Global_filter       = ifelse(.$MeanEC < 30, "fail", "pass")) 
samps$Cville_only_filter[which(samps$year == 2012 & samps$city == "Charlottesville")] = "fail"
samps$Global_filter[which(samps$year == 2012 & samps$city == "Charlottesville")] = "fail"

########################## 
# filter by consecutive year in dataset

samps %>%
  group_by(city, year) %>%
  summarise(N = n()) %>%
  dcast(city~year) -> samps_year

samps_year %<>%
  mutate(samps_per_year = rowMeans(samps_year[,-1], na.rm = T),
         NA_count= apply(., 1, function(x) sum(is.na(x))) )

samps_year %>%
  .[which(.$samps_per_year >= 1.5  &
            .$NA_count < 10),] ->
  samps_year_to_use

samps %>%
  .[which(.$city %in% samps_year_to_use$city),] %>%
 .[-which(.$year == 2012 & .$city == "Charlottesville"),]  ->
  filtered_samps_for_analysis

filtered_samps_for_analysis %>%
  .[which(.$city %in% samps_year_to_use$city),] %>%
  group_by(country) %>%
  summarise(N= n())

filtered_samps_for_analysis %>%
  .[which(.$city %in% samps_year_to_use$city),] %>%
  group_by(country, city) %>%
  summarise(N= n())

filtered_samps_for_analysis %<>% 
  mutate(Month = 
           month(as.Date(collectionDate, 
                         format = "%m/%d/%Y")))

filtered_samps_for_analysis %>% 
  ggplot(
    aes(
      y=paste(country, city, sep = " "),
      x=as.numeric(Month),
      fill = factor(month.abb[Month], levels = month.abb),
    )
  ) +
  ylab("Locale") +
  xlab("Month") +
  geom_point(size = 3.2,
             shape = 21) + 
  scale_fill_brewer(name = "Month",
                    palette = "Spectral") +
  xlim(5.5,12.5) +
  theme_bw() +
  facet_wrap(~year, nrow = 1)->
  Sample_schedule

ggsave(Sample_schedule, 
       file = "Sample_schedule.pdf",
       width = 12,
       height = 4)
ggsave(Sample_schedule, 
       file = "Sample_schedule.png",
       width = 12,
       height = 4)

#### include consecutive sample filtering
filtered_samps_for_analysis

samps %<>%
  mutate(ConsecYear_filter = ifelse(.$sampleId %in% filtered_samps_for_analysis$sampleId, "pass", "fail") ) 

samps %<>%
  mutate(ConsecCville_filter = ifelse(.$Cville_only_filter == "pass" &
                                      .$ConsecYear_filter == "pass", "pass", "fail") ) 


############################
samps_w_filters = samps #### <---------- final object for samps
############################

#Define core functions
count_mAF = function(x){
  x  %>%   
    .[which(. != 0 & . != 1)] %>%
    .[complete.cases(.)] %>%
    min ->
    mAF
  
  if(mAF == Inf){mAF = NA}
  
  return(mAF)
}

prop_missing = function(x){
  x  %>%   
    .[is.na(.)] %>%
    length() ->
    numNA
  
  return(numNA/length(x))
}

prop_fix = function(x){
  x  %>%   
    .[which(. == 0 | . == 1)] %>%
    length() ->
    Nfixed
  
  return(Nfixed/length(x))
}


# Making a global filter for all DEST populations
dat %>%
  as.data.frame() %>%
  t() %>% 
  as.data.frame -> dat_AF_t

### Make SNP filter set with Cville only
samps_w_filters %>%
  as.data.frame() %>%
  .[which(.$set == "CvilleSet"),] %>%
  .[which(.$Cville_only_filter == "pass"),] %>%
  .$sampleId ->
  VA_samps_passFilt

dat_AF_t %>%
  .[which(names(.) %in% VA_samps_passFilt )] ->
  dat_AF_t_VA

# -> describe the set
MeanAF_VA=c()
MinAF_VA=c()
Missing_VA=c()
Fixed_VA=c()

#define mean AF
apply(dat_AF_t_VA,
      1, FUN=mean, na.rm=TRUE ) -> MeanAF_VA

#define MAF counter
apply(dat_AF_t_VA,
      1, FUN=count_mAF ) -> MinAF_VA

#define and count missing sites
apply(dat_AF_t_VA,
      1, FUN=prop_missing ) -> Missing_VA

#define and count fixed sites
apply(dat_AF_t_VA,
      1, FUN=prop_fix ) -> Fixed_VA

### ### ### ### ### ### ### ### ### ### ### 
### Make SNP filter set with Global filter only
### ### ### ### ### ### ### ### ### ### ### 
samps_w_filters %>%
  as.data.frame() %>%
  .[which(.$Cville_only_filter == "pass"),] %>%
  .$sampleId ->
  Global_samps_passFilt

dat_AF_t %>%
  .[which(names(.) %in% Global_samps_passFilt )] ->
  dat_AF_t_Global

# -> describe the set
MeanAF_Global=c()
MinAF_Global=c()
Missing_Global=c()
Fixed_Global=c()

#define mean AF
apply(dat_AF_t_Global,
      1, FUN=mean, na.rm=TRUE ) -> MeanAF_Global

#define MAF counter
apply(dat_AF_t_Global,
      1, FUN=count_mAF ) -> MinAF_Global

#define and count missing sites
apply(dat_AF_t_Global,
      1, FUN=prop_missing ) -> Missing_Global

#define and count fixed sites
apply(dat_AF_t_Global,
      1, FUN=prop_fix ) -> Fixed_Global



##add to snp.dat
cbind(snp.dt, 
      MeanAF_VA, 
      MinAF_VA, 
      Missing_VA, 
      Fixed_VA,
      MeanAF_Global,
      MinAF_Global,
      Missing_Global,
      Fixed_Global
      ) ->
  snp.dt_metadata #### <---------- final object for snp.dt


############################
samps_w_filters = samps #### <---------- final object for samps
############################

## Make final object
save(samps_w_filters,
     snp.dt_metadata,
     file = "/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_Filter_Metadat.Rdata"
     )
 
