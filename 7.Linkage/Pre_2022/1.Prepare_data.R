##### ---> Correlated allele frequencies, prepping data
##### 

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)

#####

### Load data from 2L
### This is an R object congaining metadata from chr 2L
load("/project/berglandlab/jcbnunez/Shared_w_Connor/Drosophila/Cville_2L.ECfiltered.Rdata")

#### filter samples for effective coverage
filtered_samps_for_analysis %>%
  filter(MeanEC >= 30,
         city == "Charlottesville",
         year %in% 2016:2018
  ) -> samps_to_analyze

## add temperature
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"
left_join(samps_to_analyze, weather.ave) ->
  samps_metadata

#save(samps_metadata,
#     file = "samps_metadata.Rdata")

####
# Now lets create a file of allele frequencies

### open GDS file
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

### get target populations
samps <- samps_metadata

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T),
                      allele=seqGetData(genofile, "allele") )

## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId,
             snps.dt[chr%in%c("2L")]$variant.id)

snps.dt %<>% 
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))
### get allele frequency data
### 
### 
#######
dp <- seqGetData(genofile, "annotation/format/DP")
ad <- seqGetData(genofile, "annotation/format/AD")

dat <- ad$data/dp

colnames(dat) = snps.dt[chr%in%c("2L")]$SNP_id
rownames(dat) = samps$sampleId

####
###############
#Now extract the SNPs of interest
# load temperature glms
load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

glm.out %>% 
  filter(
    rnp.clean<0.01,
    mod == "aveTemp",
    chr == "2L") %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) %>%
  left_join(snps.dt) ->
  glm_2L_outP1perc

### add the inversion SNPs
inv_markers_id <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)

#### Merge both SNP ids
#### 

snp_ids_to_extract <- c(glm_2L_outP1perc$SNP_id, inv_markers_id$V1)

snp_ids_to_extract %>% length()

### extract objects $$$
dat %>%
  .[, which(colnames(dat) %in% glm_2L_outP1perc$SNP_id) ] ->
  dat_glm_p1

### Save data
###save(dat_glm_p1, glm_2L_outP1perc, samps_metadata,
###     file = "data_for_covariate_analysis_2L.Rdata")

unique_combs <- 
  combn(snp_ids_to_extract ,2) %>% 
  t() %>%
  as.data.frame()

unique_combs %>%
  separate(V1, into = c("chr1", "pos1", "feature1"), remove = F ) %>%
  separate(V2, into = c("chr2", "pos2", "feature2"), remove = F ) %>%
  mutate(bp_dist = abs(as.numeric(pos1)-as.numeric(pos2))) ->
  unique_combs_sep

save(unique_combs_sep,
     ##glm_2L_outP01perc,
     dat_glm_p1,
     file = "data_for_covariance_test.Rdata")


