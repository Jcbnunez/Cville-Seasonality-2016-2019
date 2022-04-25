#Script 10. Test the inversion model on in silico and other pools

library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
library(RColorBrewer) # customized coloring of plots
library(DescTools)   
library(rcompanion)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(tidyverse)
library(gmodels)
library(reshape2)
library(magrittr)
library(zoo)


# Load gds and metadata
ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds"
inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"

# Load prediction model
load("../SVM_2ltpred_model_and_Files.Rdata")
#model_pred_inv2lt
#inversion_markers

# Import gds
### open GDS file
genofile <- seqOpen(ingds)

### get target populations
samps <- fread(inmeta)

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

snps.dt %<>%
mutate(snp.vcf.id=paste(chr, pos, "SNP", sep = "_"))


## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

selected_snps = snps.dt %>%
                filter(snp.vcf.id %in% inversion_markers)

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=selected_snps$variant.id)

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId,
             selected_snps$variant.id)

## Extract data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")
dat <- ad$data/dp #if poolSNP
#dat <- ad/dp #if SNAPE
dim(dat)  

## Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  "SNP", sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")


dat_noNA <- na.aggregate(dat, mean)
## Predict pools by frequency
dat_noNA %>%
  as.data.frame() %>%
  rowMeans(na.rm = T) %>%
  round(2) ->
  Inversion_by_mean

#by svm
SVMpred <- predict(SVM_model_pred_inv2lt, dat_noNA) 

SVMpred %>%
  round(2) ->
  Inversion_by_SVM

data.frame(Inversion_by_mean) %>%
  cbind(Inversion_by_SVM) %>%
  as.data.frame %>%
  mutate(sampleId = rownames(.)) ->
  inversion_pred_df

######
inversion_pred_df%>%
  ggplot(
    aes(
      x=Inversion_by_mean,
      y=Inversion_by_SVM*2
    )
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  geom_smooth(method = "lm") +
  geom_point() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Inv freq. mean DEST") +
  ylab("Inv freq. SVM DEST x 2") ->
  inver_metrics

ggsave(inver_metrics,
       file = "inver_metrics.pdf")
  
#### plot
#### ### plot temp
load("./weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"

inversion_pred_df %>%
  separate(sampleId, remove = F,
           into = c("state","city","year","date_col"),
           sep ="_") %>%
  mutate(date_col_asDat = as.Date(date_col, format = "m%md%d")) %>%
  left_join(weather.ave) ->
  Inv_data_to_plot


Inv_data_to_plot %>%
  .[complete.cases(.$aveTemp),] %>%
  ggplot(aes(
    x=aveTemp,
    y=Inversion_by_mean,
    color = year
  )) +
  geom_smooth(method = "lm", se = F) +
  geom_point() +
  facet_wrap(~as.character(state)) ->
  inv_temp_plot

ggsave(inv_temp_plot,
       file = "inv_temp_plot.pdf")


Inv_data_to_plot %>%
  ggplot(aes(
    x=date_col_asDat,
    y=Inversion_by_mean,
    color = year
  )) +
  geom_line() +
  geom_point() ->
  inv_date_plot

ggsave(inv_date_plot,
       file = "inv_date_plot.pdf")

