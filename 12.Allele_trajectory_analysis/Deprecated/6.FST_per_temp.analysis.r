rm(list = ls())

#load packages
library(tidyverse)
library(magrittr)
library(forcats)
library(FactoMineR)
library(gtools)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gmodels)
library(foreach)
library(lubridate)

####
### Panel a <----------
### FST in Cville over time as a function of temperature
inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
samps <- fread(inmeta)

ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds"

####

#Generate outfile object
samps %>%
  filter(city %in% c("Charlotttesville","Charlottesville")) ->  samps_YwY_cville

samps_YwY_cville$city = gsub("Charlotttesville","Charlottesville", samps_YwY_cville$city)

print("Modify dates to as.Date format")

samps_YwY_cville$collectionDate = as.Date(samps_YwY_cville$collectionDate, 
                                      format='%m/%d/%Y')  

L = dim(samps_YwY_cville)[1]

comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

print("Create combination vector")

comp_vector %<>%
  as.data.frame() %>%
  mutate(day_diff = NA)

# Clean up names
samps_YwY_cville$city %>% unique
###
###
###



for(i in 1:dim(comp_vector)[1]) {
  
  date1=samps_YwY_cville$collectionDate[comp_vector[i,1]]
  date2=samps_YwY_cville$collectionDate[comp_vector[i,2]]
  
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  
  comp_vector$pop1[i] = samps_YwY_cville$city[comp_vector[i,1]]
  comp_vector$pop2[i] = samps_YwY_cville$city[comp_vector[i,2]]
  
  comp_vector$samp1[i] = samps_YwY_cville$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = samps_YwY_cville$sampleId[comp_vector[i,2]]
  
}

comp_vector %<>%
  mutate(city_test =
           ifelse(pop1 == pop2, "yes", "no"))



## only use observations from same population
print("Create comp vector")

comp_vector %>% 
  .[which(.$city_test == "yes"),] %>%
  #.[which(.$day_diff < 360),] %>%
  mutate(continent = 
           ifelse(.$pop1 %in% c("Charlottesville"),
                  "US", "EU" ))-> 
  comp_vector_samepops

print("End of part 3")

########################
### open GDS file
genofile <- seqOpen(ingds)

### Include DEST sets
samps <- rbind(samps[set=="CvilleSet"])


### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

#import inversion markers
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
final_in2Lt_markers %<>%
  data.frame(SNP_id = .) %>%
  separate(SNP_id, into = c("chr", "pos", "type"), sep = "_", remove = F) %>%
  mutate(win = "inv",
         class = "inv")

## Load haplotags
load("/scratch/yey2sn/Overwintering_ms/12.trajectory_analysis/haplo_tags_SNPids.Rdata")
haplo_tags_SNPids %<>% 
  filter(class == "GLM_LD") %>%
  separate(SNP_id, into = c("chr","pos", "type"), remove = F) %>%
  mutate(win = case_when(
    pos > 5155762 & pos < 5255762 ~ "win5",
    pos > 6255762 & pos < 6355762 ~ "win6",
    pos > 9505762 & pos < 9605762 ~ "win9",
    pos %in% final_in2Lt_markers$pos ~ "inv")) 

snps.dt %>% 
  filter(chr == "2L") %>%
  mutate(win = case_when(
    pos > 5155762 & pos < 5255762 ~ "win5",
    pos > 6255762 & pos < 6355762 ~ "win6",
    pos > 9505762 & pos < 9605762 ~ "win9",
    pos %in% final_in2Lt_markers$pos ~ "inv")) %>% 
  filter(pos %in% haplo_tags_SNPids$pos | win == "inv") ->
  snps.dt.haptags

#####
## choose number of alleles
snps.dt.haptags <- snps.dt.haptags[nAlleles==2]

snps.dt.haptags %<>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt.haptags$variant.id)

snps.dt.haptags[,af:=seqGetData(genofile, "annotation/info/AF")$data]

#####

### get allele frequency data
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")


###################
### Next part
### Calculate FST

  outfile = foreach(i=1:dim(comp_vector_samepops)[1], .errorhandling = "remove", .combine = "rbind")%do%{
    
    print(i/dim(comp_vector_samepops)[1] * 100)
    
    samps_to_compare = c(comp_vector_samepops$samp1[i], comp_vector_samepops$samp2[i])
    
    pool_sizes = c(samps_YwY_cville$nFlies[which(samps_YwY_cville$sampleId == comp_vector_samepops$samp1[i])],
                   samps_YwY_cville$nFlies[which(samps_YwY_cville$sampleId == comp_vector_samepops$samp2[i])])
  
  wins=c("win5", "win6", "win9", "inv")
  inner_loop =  foreach(k = 1:length(wins), .combine = "rbind" )%do%{
      
      message(wins[k])
    
      WinSNPs = filter(snps.dt.haptags, win == wins[k])$SNP_id
      
      ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ), which(colnames(ad) %in% WinSNPs ) ]
      rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ), which(colnames(dp) %in% WinSNPs ) ]
      
      pool <- new("pooldata",
                  npools=dim(ad.matrix)[1], #### Rows = Number of pools
                  nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
                  refallele.readcount=t(ad.matrix),
                  readcoverage=t(rd.matrix),
                  poolsizes=pool_sizes * 2)
      
      
      fst.out <- computeFST(pool, method = "Anova")
      
      data.frame(
        samp1=comp_vector_samepops$samp1[i],
        samp2=comp_vector_samepops$samp2[i],
        FST=fst.out$FST,
        win= wins[k]
      )
    } #window

  inner_loop

  }## i  

### import temperature
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
weather.ave -> weather.ave.1
weather.ave -> weather.ave.2

names(weather.ave.1)[1:2] = c("samp1", "aveTemp1")
names(weather.ave.2)[1:2] = c("samp2", "aveTemp2")

left_join(outfile, weather.ave.1[,1:2] ) %>% 
  left_join(weather.ave.2[,1:2]) %>% 
  mutate(temp_diff_C = abs((aveTemp1/10) - (aveTemp2/10))) %>%
  separate(samp1, into = c("state1", "city1", "year1", "date1"),  sep = "_", remove = F ) %>%
  separate(samp2, into = c("state2", "city2", "year2", "date2"),  sep = "_", remove = F  ) %>%
  mutate(year_diff = abs(as.numeric(year1)-as.numeric(year2)) ) %>%
  mutate(month_1 = month(as.POSIXlt(date1, format="m%md%d")),
         month_2 = month(as.POSIXlt(date2, format="m%md%d"))) ->
  FST_temp

save(FST_temp, file = "FST_temp.Rdata")
load("./FST_temp.Rdata")


FST_temp %>%
  filter(year_diff == 0) %>%
  ggplot(aes(
    x=temp_diff_C,
    y=FST
    #y=FST/(1-FST),
  )) +
  geom_density2d() +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  theme_bw() +
  ggtitle("Growing season FST", subtitle = "Uses only GLMs within a window") +
  facet_grid(win~.) ->
  fst_temp_fig

ggsave(fst_temp_fig, file = "fst_temp_fig.pdf", w= 3, h = 6)

#### Break up per month

FST_temp %>%
  mutate(anchor_month = case_when(month_1 != 8 ~ month_1,
                                  month_1 == 8 ~ month_1
                                  )) %>% 
  mutate(other_month = case_when(month_1 == 8 ~ month_2,
                                  month_1 != 8 ~ month_1)) %>% 
  .[complete.cases(.),] %>%
  filter(year_diff == 0) %>%
  ggplot(aes(
    x=as.factor(other_month),
    y=FST,
    color = year1
    #y=FST/(1-FST),
  )) +
  geom_boxplot() + 
  #geom_point(alpha = 0.4) +
  #geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  ggtitle("Growing season FST", subtitle = "Uses only GLMs within a window") +
  facet_grid(win~.) ->
  fst_temp_fig_month

ggsave(fst_temp_fig_month, file = "fst_temp_fig_month.pdf", w= 9, h = 6)




