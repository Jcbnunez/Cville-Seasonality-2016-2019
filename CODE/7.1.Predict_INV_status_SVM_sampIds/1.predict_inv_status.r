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
    SVM_pred <= 0.35 ~ "STD",
    SVM_pred >= 0.75 ~ "INV",
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

INV_2lt_markers_w_Allele = colnames(training_data)[-95]

save(sample_metadata,
     INV_2lt_markers_w_Allele,
     SVM_model_pred_inv2lt.new,
     file = "/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM.toPredict_INV_STD.new.Apr1.2022.Rdata"
     )

#load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM.toPredict_INV_STD.new.Apr1.2022.Rdata")

sample_metadata %>%
  filter(INV_STATUS %in% c("INV", "STD")) %>%
  group_by(pop,INV_STATUS ) %>%
  summarise(N= n()) %>%
  dcast(pop~INV_STATUS) %>%
  as.data.frame()

### Write sample file name
sample_metadata %>%
  filter(pop != "NY") %>%
  .$sampleId -> samples_hom_in2lt



write.table(samples_hom_in2lt, 
            file = "/project/berglandlab/Dmel_Single_Individuals/in2Lt_status_info/samples_hom_in2lt.ids.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


#####