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
library(vcfR)
library(forcats)
library(FactoMineR)
library(adegenet)
library(Hmisc)

### GENERATE POOL OF RANMLY SAMPLED SNPS in DGRP
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")

##load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/DGRP_2lt_Markers##/imputated_data_DGRP_dat.2L.Rdata")
##set.seed(123456)
##imputated_data_DGRP %>% 
##  rownames() %>%
##  .[grep("SNP", .)] %>% 
##  .[sample(length(.), 5000, replace = F)] ->
##  sampled_SNPs
##
####include inversion markers
##write.table(unique(c(sampled_SNPs, final_in2Lt_markers)), 
##            file = "sample_plus_invMarks_SNPs.2l.txt", 
##            append = FALSE, quote = FALSE, sep = "\t",
##            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
##            col.names = FALSE, qmethod = c("escape", "double"),
##            fileEncoding = "")
##
### After having created VCF subsets, procceed here
## load datasets
setwd("/scratch/yey2sn/Overwintering_ms/15.ME_PA_extra_pops")
invst <- fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
names(invst)[1] = "DGRP_Line"
invst$line_name = gsub("DGRP" , "line" , invst$DGRP_Line)

## load the DGRP
DGRP_invM <- read.vcfR(
  "./DGRP.2L.Rand_INV_Markers.recode.vcf.gz")
DGRP_invMgl <- vcfR2genind(DGRP_invM, return.alleles = T) 
DGRP_invMgl_tb = tab(DGRP_invMgl, NA.method = "asis")



###
### Train Model
DGRP_invMgl_tb %>%
  as.data.frame() %>% 
  .[complete.cases(.),] %>% 
  .[-which(rownames(.) %in% c("line_492", "line_48") ),] %>% 
  mutate(line_name = rownames(.)) %>% 
  left_join(invst[,c("line_name","In_2L_t")]) %>%
  mutate(class = case_when(
    In_2L_t == "ST" ~ 0,
    In_2L_t == "INV" ~ 1
  )) %>% 
  select(!c(line_name, In_2L_t)) -> training_data

tune.out <- tune(svm,
                 class ~., 
                 data = training_data, 
                 kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

SVM_model_pred_inv2lt.new <- tune.out$best.model

#### Cville

CVILLE_invM <- read.vcfR(
  "./CM.2L.Rand_INV_Markers.recode.vcf.gz")
CVILLE_invMgl <- vcfR2genind(CVILLE_invM, return.alleles = T) 
CVILLE_invMgl_tb = tab(CVILLE_invMgl, NA.method = "asis")


CVille_predict = CVILLE_invMgl_tb[,colnames(training_data)[-95]]

SVMpred_pa_ln.CVILLE <- predict(SVM_model_pred_inv2lt.new, 
                         CVille_predict)

data.frame(SVM_pred = SVMpred_pa_ln.CVILLE,
           sampleId = names(SVMpred_pa_ln.CVILLE)) ->
  SVMpred_pa_ln.Cville
#### Plot histogram
SVMpred_pa_ln.Cville %>%
  ggplot(aes(SVM_pred)) +
  geom_histogram() ->
  SVM_histogram
ggsave(SVM_histogram, file = "SVM_histogram.Cville.pdf")

### MAKE PCA
CVille_predict %>%
  PCA(graph = F, ncp = 5) -> CVille_predict_pca

##plot pca
CVille_predict_pca$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(SVMpred_pa_ln.Cville) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = SVM_pred
  )) +
  geom_jitter(
    width = 0.1,
    height = 0.1,
    shape = 21,
    size = 4
  ) +
  scale_fill_gradient2(low="red", high="blue", midpoint = 0.5) +
  ggtitle("CVille Samps") ->
  Cville.PCA_chr2_snps

ggsave(Cville.PCA_chr2_snps, file = "Cville.PCA_chr2_snps.pdf")
####

####
ME_PA_invM <- read.vcfR(
  "./Taylors.2L.Rand_INV_Markers.recode.vcf.gz")
ME_PA_invMgl <- vcfR2genind(ME_PA_invM, return.alleles = T) 
ME_PA_invMgl_tb = tab(ME_PA_invMgl, NA.method = "asis", )

ME_PA_predict = ME_PA_invMgl_tb[,colnames(training_data)[-95]]

all_column_mean <- apply(ME_PA_predict, 2, mean, na.rm=TRUE)
for(k in 1:length(colnames(ME_PA_predict)) ){
  
  i=colnames(ME_PA_predict)[k]
  ME_PA_predict[,i][is.na(ME_PA_predict)[,i]] <- all_column_mean[i]
  
}


SVMpred_pa_ln.WORLD <- predict(SVM_model_pred_inv2lt.new, 
                               ME_PA_predict)

data.frame(SVM_pred = SVMpred_pa_ln.WORLD,
           sampleId = names(SVMpred_pa_ln.WORLD)) ->
  SVMpred_pa_ln.WORLD

SVMpred_pa_ln.WORLD %>%
  ggplot(aes(SVM_pred)) +
  geom_histogram() ->
  SVM_histogram.World
ggsave(SVM_histogram.World, file = "SVM_histogram.Wolrd.pdf")

### MAKE PCA
ME_PA_predict %>%
  PCA(graph = F, ncp = 5) -> ME_PA_predict_pca

##plot pca
ME_PA_predict_pca$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(SVMpred_pa_ln.WORLD) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = SVM_pred
  )) +
  geom_jitter(
    width = 0.1,
    height = 0.1,
    shape = 21,
    size = 4
  ) +
  scale_fill_gradient2(low="red", high="blue", midpoint = 0.5) +
  ggtitle("Taylor Samps") ->
  Wolrd.PCA_chr2_snps

ggsave(Wolrd.PCA_chr2_snps, file = "Wolrd.PCA_chr2_snps.pdf")

#####
samples_set = c(rownames(SVMpred_pa_ln.Cville), rownames(SVMpred_pa_ln.WORLD))

rbind(
SVMpred_pa_ln.Cville,
SVMpred_pa_ln.WORLD)  %>% 
  mutate(INV_STATUS = case_when(
    SVM_pred <= 0.75 ~ "STD",
    SVM_pred >= 0.35 ~ "INV",
    SVM_pred > 0.35 & SVM_pred < 0.75 ~ "HET"
  )) %>%
  mutate(pop = case_when(
    sampleId %in% samples_set[grep("CM", samples_set)] ~ "CM",
    sampleId %in% samples_set[grep("OW", samples_set)] ~ "CM",
    sampleId %in% samples_set[grep("line", samples_set)] ~ "DGRP",
    sampleId %in% samples_set[grep("ME", samples_set)] ~ "ME",
    sampleId %in% samples_set[grep("LN", samples_set)] ~ "PA",
    sampleId %in% samples_set[grep("CO", samples_set)] ~ "Cameroon",
    sampleId %in% samples_set[grep("FR", samples_set)] ~ "France",
    sampleId %in% samples_set[grep("GU", samples_set)] ~ "Guinea",
    sampleId %in% samples_set[grep("^I", samples_set)] ~ "NY",
    sampleId %in% samples_set[grep("Netherlands", samples_set)] ~ "Netherlands",
    sampleId %in% samples_set[grep("ZI", samples_set)] ~ "Zambia"
  )) -> sample_metadata
sample_metadata$pop[is.na(sample_metadata$pop)] = "Carib" 

sample_metadata %>%
  filter(INV_STATUS %in% c("INV", "STD")) %>%
  group_by(pop,INV_STATUS ) %>%
  summarise(N= n()) %>%
  as.data.frame()

#tab(DGRP_invMgl)[,grep("\\.0", colnames(tab(DGRP_invMgl)))] -> DGRP_invMgl_tb
#colnames(DGRP_invMgl_tb) = gsub("\\.0", "", colnames(DGRP_invMgl_tb))


#tab(CVILLE_invMgl)[,grep("\\.0", colnames(tab(CVILLE_invMgl)))] -> CVILLE_invMgl_tb
#colnames(CVILLE_invMgl_tb) = gsub("\\.0", "", colnames(CVILLE_invMgl_tb))

intersect(intersect(colnames(CVILLE_invMgl_tb),
                    colnames(DGRP_invMgl_tb)),
          colnames(ME_PA_invMgl_tb)) -> sampled_SNPs

### Make joint object
rbind(
CVILLE_invMgl_tb[,sampled_SNPs],
DGRP_invMgl_tb[,sampled_SNPs],
ME_PA_invMgl_tb[,sampled_SNPs]) ->
  joint_sample_2l_snps

invst$line_name[-which(invst$In_2L_t %in% c("ST","INV") )] -> remove_DGRP_hets

data.frame(sampleId = rownames(joint_sample_2l_snps))

final_in2Lt_markers[final_in2Lt_markers %in% colnames(joint_sample_2l_snps)] ->
  shared_inv_markers



####
predict_data =  joint_sample_2l_snps[,shared_inv_markers] 

### impute missing data
all_column_mean <- apply(predict_data, 2, mean, na.rm=TRUE)
for(k in 1:length(colnames(predict_data)) ){
  
  i=colnames(predict_data)[k]
  predict_data[,i][is.na(predict_data)[,i]] <- all_column_mean[i]
  
}
#### predict
SVMpred_pa_ln <- predict(SVM_model_pred_inv2lt, 
                         predict_data
) 
data.frame(SVM_pred = SVMpred_pa_ln,
           sampleId = names(SVMpred_pa_ln)) ->
  SVMpred_pa_ln
#### Plot histogram
SVMpred_pa_ln %>%
  ggplot(aes(SVM_pred)) +
  geom_histogram() ->
  SVM_histogram
ggsave(SVM_histogram, file = "SVM_histogram.pdf")

## add labels
SVMpred_pa_ln %<>%
mutate(INV_STATUS = case_when(
  SVM_pred >= 0.70 ~ "STD",
  SVM_pred <= 0.30 ~ "INV",
  SVM_pred > 0.30 & SVM_pred < 0.70 ~ "HET"
))


### MAKE PCA
predict_data %>%
  PCA(graph = F, ncp = 5) -> predict_data_pca

##plot pca
predict_data_pca$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(., sample_metadata ) %>%
  left_join(SVMpred_pa_ln) -> pca_analysis_metadat
######
pca_analysis_metadat %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = SVM_pred
  )) +
  geom_jitter(
    width = 0.1,
    height = 0.1,
    shape = 21,
    size = 4
  ) +
  scale_fill_gradient2(low="red", high="blue", midpoint = 0.5) +
  facet_wrap(~pop) +
  ggtitle("PCA with Inversion informative Markers ONLY") ->
  PCA_chr2_snps

ggsave(PCA_chr2_snps, file = "PCA_chr2_snps.pdf")
####
####

###### Create the list of inversion or derived
pca_analysis_metadat %<>%
  mutate(INV_STATUS = case_when(
    SVM_pred >= 0.70 ~ "STD",
    SVM_pred <= 0.30 ~ "INV",
    SVM_pred > 0.30 & SVM_pred < 0.70 ~ "HET"
  ))
####
pca_analysis_metadat %>%
  filter(INV_STATUS %in% c("INV", "STD")) %>%
  group_by(pop,INV_STATUS ) %>%
  summarise(N= n()) %>%
  as.data.frame()




#### load prediction model
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")



tab(ME_PA_invMgl)[,grep("\\.0", colnames(tab(ME_PA_invMgl)))] -> ME_PA_invMgl_tb
colnames(ME_PA_invMgl_tb) = gsub("\\.0", "", colnames(ME_PA_invMgl_tb))

data.frame(marker=final_in2Lt_markers,
           TF=final_in2Lt_markers %in% colnames(ME_PA_invMgl_tb) 
           )
### Impute data
# getting median of each column using apply() 
all_column_median <- apply(as.data.frame(ME_PA_invMgl_tb), 2, mean, na.rm=TRUE)
data_f=as.data.frame(ME_PA_invMgl_tb)
  
# imputing median value with NA 

SVMpred_pa_ln %>%
  #round(1) %>%
  data.frame(SVM_pred = .) %>%
  mutate(sampleId= rownames(.)) ->
  SVMpred_pa_ln_svm

#make PCA
ME_PA_invMgl_tb_noNA %>%
  PCA(graph = F, ncp = 5) -> ME_PA_invMgl_tb_noNA_pca

SVMpred_pa_ln_pca$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId= rownames(.)) %>%
  left_join(SVMpred_pa_ln_svm) -> prediction_taylors_dat


prediction_taylors_dat %>%
  mutate(pop = case_when(sampleId %in% prediction_taylors_dat$sampleId[grep("ME", prediction_taylors_dat$sampleId )] ~ " ME",
                         sampleId %in% prediction_taylors_dat$sampleId[grep("LN", prediction_taylors_dat$sampleId )] ~ " LN" )           
         ) -> prediction_taylors_dat

prediction_taylors_dat %>%
  ggplot(
    aes(
      SVM_pred,
      fill=pop
    )
  ) +
  geom_histogram() ->
  SVM_hist
ggsave(SVM_hist,
       file = "SVM_hist.pdf")



#####
#####
#####
#####