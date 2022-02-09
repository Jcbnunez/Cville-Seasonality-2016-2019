### This script prepares data for haplovalidate

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)

### Load data from 2L
load("/project/berglandlab/jcbnunez/Shared_w_Connor/Drosophila/Cville_2L.ECfiltered.Rdata")

#### filter samples for effective coverage
filtered_samps_for_analysis %>%
  filter(MeanEC >= 30,
         city == "Charlottesville",
         year %in% 2016:2018
         ) -> samples_for_haplovalidate

#write.table(samples_for_haplovalidate, 
#            file = "sampled_for_haplovalidate.all.txt", 
#            append = FALSE, 
#            quote = FALSE, 
#            sep = "\t",
#            eol = "\n", 
#            na = "NA", 
#            dec = ".", 
#            row.names = FALSE,
#            col.names = TRUE, qmethod = c("escape", "double"),
#            fileEncoding = "")
#
#### start preparing the sync file


##################################
# Part 1 -- Prepare Files
##################################

### open GDS file
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

### get target populations
samps <- samples_for_haplovalidate

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

### get allele frequency data
### 
### 
#######
dp <- seqGetData(genofile, "annotation/format/DP")
dp_t <- t(dp)
colnames(dp_t) = samps$sampleId
dp_t <- cbind(snps.dt[,c(1:2,6)], dp_t)

dp_t %>%
  separate(allele, 
           into = c("REF", "ALT"),
           sep = ",") %>% 
  melt(id = c("chr","pos","REF","ALT"), value.name = "Coverage") ->
  dp_melt

########
ad <- seqGetData(genofile, "annotation/format/AD")
ad_t <- t(ad$data)
colnames(ad_t) = samps$sampleId
ad_t <- cbind(snps.dt[,c(1:2,6)], ad_t)

ad_t %>%
  separate(allele, 
           into = c("REF", "ALT"),
           sep = ",") %>% 
  melt(id = c("chr","pos","REF","ALT"), value.name = "ALT_count") ->
  ad_melt

####
#######
cbind(dp_melt, ALT_count = ad_melt$ALT_count) %>% 
  mutate(REF_count = Coverage-ALT_count) -> ad_dp_obj

save(ad_dp_obj,
     file = "ad_dp_obj.Rdata")

ad_dp_obj %>%
dcast(chr+pos+variable~REF, value.var = "REF_count" ) -> ref_count
ref_count[is.na(ref_count)] <- 0
names(ref_count)[4:7] = gsub("$","_ref", names(ref_count)[4:7])

ad_dp_obj %>%
dcast(chr+pos+variable~ALT, value.var = "ALT_count" ) -> alt_count
alt_count[is.na(alt_count)] <- 0
names(alt_count)[4:7] = gsub("$","_alt", names(alt_count)[4:7])

#####
#####

cbind(ref_count, alt_count[,4:7]) %>% 
  mutate(A = A_ref+A_alt,
         T = T_ref+T_alt,
         C = C_ref+C_alt,
         G = G_ref+G_alt,
         N = 0,
         DEL = 0
         ) %>%
  mutate(sync_form = paste(A,T,C,G,N,DEL, sep = ":")) -> sync_format_table

save(sync_format_table,
     file = "sync_format_table_interm_file.Rdata")


