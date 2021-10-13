## This code is used to delineate the 2Lt inversion using the DGRP
## 

#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)

#Metadata
invst <- fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
names(invst)[1] = "DGRP_Line"
invst$line_name = gsub("DGRP" , "line" , invst$DGRP_Line)

# Load the DGRP vcf
DGRP <- read.vcfR(
  "/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz")

#transform to genlight object
DGRPgl <- vcfR2genlight(DGRP) 
#check individuals
#tab(DGRPgl) %>%
#  rownames()

#Extract data for imputation
tab(DGRPgl, NA.method = "asis") %>%
  as.data.frame() %>% 
  t() ->
  all_dgrp_data
  
#Separate data based on inversion set
standard  <- invst$line_name[which(invst$In_2L_t == "ST")]
inversion <- invst$line_name[which(invst$In_2L_t == "INV")]
heterokaryotype <- invst$line_name[which(invst$In_2L_t == "INV/ST")]

#Make standard set
all_dgrp_data %>%
  .[,which(colnames(.) %in% standard)] ->
  standard_lines

all_dgrp_data %>%
  .[,which(colnames(.) %in% inversion)] ->
  inverted_lines

all_dgrp_data %>%
  .[,which(colnames(.) %in% heterokaryotype)] ->
  heterokaryotype_lines

######## Define functions to study sets

filter_drgp_sets = function(x, miss_thresh){
  
  count_NA = function(x) {
    x  %>%   
      .[is.na(.)] %>%
      length() ->
      numNA
        return(numNA/length(x))
  }
  
  apply(x,
        2, 
        FUN=count_NA
  ) -> x_NA_lines
  
  x_keep_samps <- names(x_NA_lines[x_NA_lines <= miss_thresh])
  
  x %>%
    .[,which(colnames(.) %in%  x_keep_samps)] -> 
    x_lines_samps_filt
  
  #count NA per SNPs
  apply(x_lines_samps_filt,
        1, 
        FUN=count_NA
  ) -> x_NA_SNPs
  
  x_keep_snps <- names(x_NA_SNPs[x_NA_SNPs <= miss_thresh])
  
  x_lines_samps_filt %>%
    .[which(rownames(.) %in%  x_keep_snps),] -> 
    x_lines_samps_SNPs_filt
  
  x_lines_samps_SNPs_filt_NAimp = na.aggregate(x_lines_samps_SNPs_filt,
                                               FUN= mean)
  
  return(x_lines_samps_SNPs_filt_NAimp)
  
}

######## ######## ######## ######## ######## 
######## Evaluate the standard_lines
standard_lines_flt = filter_drgp_sets(standard_lines, 0.05)

inverted_lines_flt = filter_drgp_sets(inverted_lines, 0.05)

heterokaryotype_lines_flt = filter_drgp_sets(heterokaryotype_lines, 0.25)

### Merge the datasets
standard_lines_flt %<>% 
  as.data.frame() 
  names(standard_lines_flt) = gsub("^", "STD_", names(standard_lines_flt))
  standard_lines_flt %<>% mutate(SNP_id = rownames(.))
  
inverted_lines_flt %<>% 
  as.data.frame()
  names(inverted_lines_flt) = gsub("^", "INV_", names(inverted_lines_flt))
  inverted_lines_flt %<>% mutate(SNP_id = rownames(.))
  
heterokaryotype_lines_flt %<>% 
  as.data.frame() 
  names(heterokaryotype_lines_flt) = gsub("^", "HET_", names(heterokaryotype_lines_flt))
  heterokaryotype_lines_flt %<>% mutate(SNP_id = rownames(.))
  
right_join(standard_lines_flt, inverted_lines_flt, by = "SNP_id") %>% 
  right_join(heterokaryotype_lines_flt, by = "SNP_id" ) ->
  imputated_data_DGRP

SNP_id_col = which(names(imputated_data_DGRP) == "SNP_id")

imputated_data_DGRP %>%
  .[grep("SNP" ,imputated_data_DGRP$SNP_id),] %>% 
  .[,-SNP_id_col] ->
  imputated_data_DGRP_dat

rownames(imputated_data_DGRP_dat) = imputated_data_DGRP[,SNP_id_col][grep("SNP" ,imputated_data_DGRP$SNP_id)]


#Save object
save(imputated_data_DGRP_dat,
     file = "./imputated_data_DGRP_dat.2L.Rdata")

########
######## Now build PCA
########



imputated_data_DGRP_dat %>% 
  t() %>%
  as.data.frame() %>% 
  PCA(scale.unit = FALSE, 
      ncp = 5,
      graph = F
      #ind.sup = grep("HET", names(imputated_data_DGRP_dat))
      ) -> 
  PCA_obj_2LT

#save(PCA_obj_2LT, file = "PCA_obj_2LT.Rdata")
#load("./PCA_obj_2LT.Rdata")

PCA_obj_2LT$ind$coord %>% 
  as.data.frame() %>%
  mutate(DGRP_Line = rownames(.)) ->
  PCA_projs

#cluster analysis -- preeliminary
# K-Means Cluster Analysis
set.seed(1234)
fit <- kmeans(PCA_projs[,c(1,2)], 3) # 5 cluster solution
# append cluster assignment
PCA_projs <- data.frame(PCA_projs, fit$cluster)
PCA_projs$DGRP_Line_Name = PCA_projs$DGRP_Line
PCA_projs$DGRP_Line = gsub("line", "DGRP", PCA_projs$DGRP_Line)
PCA_projs %>%
  left_join(invst) ->
    PCA_projs
  
  PCA_projs %>%
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.2,
      color = perc_comp,
      shape = as.factor(fit.cluster))
  ) +
  geom_point(size = 3) -> 
  PCA_graph_2lt
  
ggsave(PCA_graph_2lt,
       file ="PCA_graph_2lt.pdf")

## Remove data not used to train
## Concordance test

PCA_projs$test = NA

PCA_projs$test[PCA_projs$fit.cluster == 3] = "fail"
PCA_projs$test[PCA_projs$fit.cluster == 1 & PCA_projs$In_2L_t == "ST"] = "pass"
PCA_projs$test[PCA_projs$fit.cluster == 1 & PCA_projs$In_2L_t != "ST"] = "fail"
PCA_projs$test[PCA_projs$fit.cluster == 2 & PCA_projs$In_2L_t == "INV"] = "pass"
PCA_projs$test[PCA_projs$fit.cluster == 2 & PCA_projs$In_2L_t != "INV"] = "fail"

keep_samps <- PCA_projs$DGRP_Line_Name[which(PCA_projs$test == "pass")]

keep_samps_homoKaritypes = keep_samps
save(keep_samps_homoKaritypes,
     file = "./DGRP_homoKariotypes.Rdata")

### Repat PCA
#Make PC
tab(DGRPgl) %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% keep_samps),] %>%
  PCA(scale.unit = FALSE, 
      ncp = 5,
      graph = F) -> 
  PCA_obj_2LT_filt


PCA_obj_2LT_filt$ind$coord %>% 
  as.data.frame() %>%
  mutate(DGRP_Line = gsub("line","DGRP",rownames(.)) ,
         DGRP_Line_Name = rownames(.)) %>%
  left_join(invst) ->
  PCA_projs_filt


PCA_projs_filt %>%
  ggplot(
    aes(
      x=Dim.1,
      y=Dim.2,
      color = In_2L_t)
  ) +
  geom_point(size = 3) -> 
  PCA_graph_2lt_flt

ggsave(PCA_graph_2lt_flt,
       file ="PCA_graph_2lt_flt.pdf")

## Save object:
save(PCA_obj_2LT_filt,
     file = "PCA_obj_2LT_filt.Rdata")

### Run correlation analysis
### MUST RUN VIA JOB ... takes ages otherwise!
##dimdesc(PCA_obj_2LT_filt,
##        axes = 1:3, 
##        proba = 1
##        ) -> PCA_graph_2lt_flt_dimdesc

load("./PCA_graph_2lt_flt_dimdesc.Rdata")

PCA_graph_2lt_flt_dimdesc$Dim.1 %>% 
  as.data.frame() %>%
  .[complete.cases(.),] %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, into = c("chr","pos","Type"), remove = F) %>%
  mutate(quanti.p.value.adj = p.adjust(quanti.p.value, method = "bonferroni")) %>% 
  mutate(Cor_value = ifelse(.$quanti.correlation^2 > 0.9, "informative", "Cor<NotInformative%" )) ->
    correlation_info

correlation_info$chr = gsub("X","",correlation_info$chr)  
correlation_info$SNP_id = gsub("X","",correlation_info$SNP_id)  

correlation_info %>%
  group_by(Cor_value) %>%
  summarise(N = n())
# A tibble: 2 × 2
# P_thresh1      N
# 1 Cor<70%   998028
# 2 Cor>70%     8040

#save object
save(correlation_info,
     file = "./Inv2lt_informative_markers.Rdata")

##########
#Make SNP vector
selected_SNPs <- correlation_info$SNP_id[which(correlation_info$Cor_value == "informative")]
##########


# Define in2Lt traditional boundaries
in2lt_beg=2225744	
in2lt_end=13154180


correlation_info %>%
  ggplot(
    aes(
      x=as.numeric(pos),
      y=-log10(quanti.p.value.adj),
      color = Cor_value
    )
  ) + geom_point() +
  geom_vline(xintercept = in2lt_beg) +
  geom_vline(xintercept = in2lt_end) ->
  corr_plot_2lt

ggsave(corr_plot_2lt,
       file = "corr_plot_2lt.png",
       width = 6,
       height = 3)



###########
############ Train DAPC model  
load("./Inv2lt_informative_markers.Rdata")
selected_SNPs <- correlation_info$SNP_id[which(correlation_info$Cor_value == "informative")]
selected_SNPs  ## <-- comes from above 
load("./DGRP_homoKariotypes.Rdata")
keep_samps_homoKaritypes

#Make new PC with informative SNPs and samples
#First make a SNP object with the training samples
tab(DGRPgl) %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% keep_samps_homoKaritypes),] %>% 
  t() %>%
  as.data.frame() %>% 
  .[which(rownames(.) %in% selected_SNPs),] %>%
  t() %>%
  as.data.frame() -> 
  Inv_inf_SNPs_df

#Make 
prop_missing = function(x){
  x  %>%   
    .[which(. %in% c(0,1,2) )] %>%
    length() ->
    numNA
  
  return(numNA/length(x))
}


#apply function
Inv_inf_SNPs_df%>%
  as.data.frame() %>%
  apply(.,
        1, FUN=prop_missing ) -> 
  prop_missing_homoka

Inv_inf_SNPs_df %<>%
  mutate(line_name = rownames(.))

left_join(invst, 
          data.frame(line_name = names(prop_missing_homoka) , homoka_com = prop_missing_homoka)
          ) %>%
left_join(Inv_inf_SNPs_df) %>%
  melt(id = c("DGRP_Line",
      "In_2L_t",
      "In_2R_NS",
      "In_2R_Y1",
      "In_2R_Y2",
      "In_2R_Y3",
      "In_2R_Y4",
      "In_2R_Y5",
      "In_2R_Y6",
      "In_2R_Y7",
      "In_3L_P",
      "In_3L_M",
      "In_3L_Y",
      "In_3R_P",
      "In_3R_K",
      "In_3R_Mo",
      "In_3R_C",
      "perc_comp",
      "line_name",
      "homoka_com"
      )) %>% 
  separate(variable,
           into = c("chr", "pos", "type"), sep = "_") ->
  inversion_markers_guides
  
  #filter missing data
  inversion_markers_guides %>%
    .[,c("In_2L_t",
        "perc_comp",
        "line_name",
        "homoka_com",
        "chr",
        "pos",
        "type",
        "value")] %>% 
    .[which(.$value %in% c(0,1,2) ),] %>% 
    dcast(line_name+In_2L_t~pos, value.var = "value") %>% 
    .[complete.cases(.),] %>% .$In_2L_t



  ggplot(
    aes(
      x=as.numeric(pos),
      y=value,
      color = In_2L_t
    )
  ) + 
  geom_point() ->
  inversion_markers


ggsave(inversion_markers,
       file = "inversion_markers.pdf")




Inv_inf_SNPs_df %>% dim()
#Second make a SNP object with the test samples. Which are all heterozygiotes
invst$In_2L_t %>% table
invst$DGRP_Line[which(invst$In_2L_t == "INV/ST")] -> het_samples
het_samples = gsub("DGRP","line", het_samples)

tab(DGRPgl) %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% het_samples),] %>% 
  t() %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% selected_SNPs),] %>%
  t() %>%
  as.data.frame() -> 
  Inv_inf_SNPs_Heterozygotes

Inv_inf_SNPs_Heterozygotes %>% dim()

### Train DAPC model
data.frame(DGRP_Line = gsub("line", "DGRP", rownames(Inv_inf_SNPs_df))) %>%
  left_join(invst) -> sample_indentities

grp <- as.factor(sample_indentities$In_2L_t)

#png("crossval.png")
#xval <- xvalDapc(Inv_inf_SNPs_df, grp, 
#                 n.pca.max = 300, training.set = 0.9,
#                 result = "groupMean", center = TRUE, scale = FALSE,
#                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
#dev.off()

dapc_2lt <- dapc(Inv_inf_SNPs_df, grp,
                 n.pca=1,n.da=1)

#probability of agreement
pred.sup <- predict.dapc(dapc_2lt, newdata=Inv_inf_SNPs_Heterozygotes)





