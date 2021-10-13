#Load libraries
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
library(adegenet)
library(vcfR)
library(forcats)
library(FactoMineR)


#Metadata
invst <- fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
names(invst)[1] = "DGRP_Line"
invst$line_name = gsub("DGRP" , "line" , invst$DGRP_Line)

# Load the DGRP vcf
DGRP_invM <- read.vcfR(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/inv2Lt_markers.DGRP.2L.recode.vcf.gz")
DGRPinvgl <- vcfR2genlight(DGRP_invM) 

CM_invM <- read.vcfR(
  "/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/inv2Lt_markers.CM.2L.recode.vcf.gz")
CMinvgl <- vcfR2genlight(CM_invM) 

#Extract data for imputation
tab(DGRPinvgl, NA.method = "asis") %>%
  as.data.frame() %>%
  .[complete.cases(.),] %>%
  t() ->
  DGRPinvgl_dat

tab(CMinvgl, NA.method = "asis") %>%
  as.data.frame() %>%
  .[complete.cases(.),] %>%
  t() ->
  CMinvgl_dat

tab(CMinvgl, NA.method = "asis") %>%
  as.data.frame() %>%
  .[complete.cases(.),]  ->
  CMinvgl_dat_t

####
DGRPinvgl_dat %>%
  .[,-which(colnames(.) %in% c("line_492", "line_48") )] ->
  DGRPinvgl_dat_flt


  DGRPinvgl_dat_flt %>%
  colnames(.) %>%
  data.frame(line_name = .) %>%
  left_join(invst)  %>% 
  select(line_name, In_2L_t) ->
  dapc_metadat

  
  ## data for training
  ## 
  DGRPinvgl_dat_flt %>%
    t() %>%
    as.data.frame() ->
    DGRPinvgl_dat_flt_t
  
  
  dapc_metadat$In_2L_t %>% table
  #.
  #INV  ST 
  #19 121 
  #DGRPinvgl_dat_flt_t_H = rbind(DGRPinvgl_dat_flt_t, inSilicoHet= rep(1, 47)) 
  
  inversion_markers <-  rownames(DGRPinvgl_dat_flt) #3 <<< --- inv markers

#### Determine inversion markers in the DEST data
  
  
  ####### Validate frequencies in DEST
  ####### 
  
  # Load gds and metadata
  ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds"
  inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
  
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
  
  #########
  dat %>%
    as.data.frame() %>%
    mutate( sampleId = rownames(.)) %>% 
    melt(id = "sampleId") %>% 
    left_join(samps) %>% 
    group_by(variable, set) %>%
    summarize(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm = T)) %>%
ggplot(aes(
  x=fct_reorder(variable, mean),
  y=mean,
  ymin=mean-sd,
  ymax=mean+sd
)) +
    geom_errorbar() +
    geom_point() +
    coord_flip() +
    facet_wrap(~set) +
    theme(axis.text.y = element_text(size = 6)) -> 
    marker_frequencies.set
  
  ggsave(marker_frequencies.set,
         file = "marker_frequencies.set.png")
    
  
  dat %>%
    as.data.frame() %>%
    mutate( sampleId = rownames(.)) %>% 
    melt(id = "sampleId") %>% 
    left_join(samps) %>% 
    group_by(variable, country) %>%
    summarize(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm = T)) %>%
    ggplot(aes(
      x=fct_reorder(variable, mean),
      y=mean,
      ymin=mean-sd,
      ymax=mean+sd
    )) +
    geom_errorbar() +
    geom_point() +
    coord_flip() +
    facet_wrap(~country) +
    theme(axis.text.y = element_text(size = 6)) -> 
    marker_frequencies.country 
  
  ggsave(marker_frequencies.country,
         file = "marker_frequencies.country.png",
         width = 8,
         height = 8)
  
### Generate final list of markers  
####
  dat %>% colnames -> final_in2Lt_markers
####
###  
  
  train_data <- DGRPinvgl_dat_flt_t[,which(names(DGRPinvgl_dat_flt_t) %in%  final_in2Lt_markers)] 
    
  train_data$class = ifelse(dapc_metadat$In_2L_t == "ST", 0, 1)

  tune.out <- tune(svm,
                   class ~., 
                   data = train_data, 
                   kernel = "linear",
                   ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

  # extract the best model
  SVM_model_pred_inv2lt <- tune.out$best.model
  
  ypred <- predict(SVM_model_pred_inv2lt, CMinvgl_dat_t)
  
  data.frame(SVM = round(ypred,2), line_name = names(ypred)) ->
    CM_predictions
  
  save(SVM_model_pred_inv2lt,
       CM_predictions,
       final_in2Lt_markers,
       file = "SVM_2ltpred_model_and_Files.Rdata")
  

  #### plot
  ####
  rbind(DGRPinvgl_dat_flt_t[,which(names(DGRPinvgl_dat_flt_t) %in%  final_in2Lt_markers)],
        CMinvgl_dat_t[,which(names(CMinvgl_dat_t) %in%  final_in2Lt_markers)]) -> combined_dat 
  
  ###  
  combined_dat %>%
    t() %>%
    PCA(graph = F, scale.unit = F) ->
    dgrp_mark_pca
  
  dgrp_mark_pca$ind$coord %>%
    as.data.frame() %>%
    mutate(line_name = rownames(.)) %>%
    left_join(invst[,c( "line_name", "In_2L_t" )]) %>%
  left_join(CM_predictions) %>%
  mutate(SVM= ifelse(is.na(.$In_2L_t), SVM, ifelse(.$In_2L_t == "ST", 0, 1)  )) -> pca_annot
  
  pca_annot[which(pca_annot$In_2L_t %in% c("INV", "ST")),"group"] = "DGRP"
  pca_annot[is.na(pca_annot$In_2L_t), "group"] = "Unknown (C. Mt.)"
  
  pca_annot %>%
    ggplot(aes(
      x=Dim.1,
      y=Dim.2,
      fill=SVM
    )) +
    geom_jitter(size = 3,
               alpha = 0.9,
               shape = 21,
               width = 0.05,
               height = 0.05)+
    ggtitle("Support Vector Machine Model") +
    scale_fill_gradient2(name= "Inv. Classification",
                         low = "blue", 
                        high = "red", 
                        midpoint = 0.5, 
                        na.value = NA) +
    facet_wrap(~group, ncol = 2) +
    theme(legend.position =  "bottom") +
    theme_bw() ->
    dgrp_2ltm_graph
  
  ggsave(dgrp_2ltm_graph, 
         file = "dgrp_2ltm_graph.pdf",
         width = 6,
         height = 3)
  
  ### table
  pca_annot %>%
    filter(group != "DGRP") %>%
    separate(line_name,
             into = c("exp", "id", "date"),
             remove = F, sep = "_") %>%
    filter(exp %in% c("OW1","OW2")) %>%
    mutate(SVM_t = round(SVM,1)) %>%
    filter(SVM_t %in% c(0,0.5,1) ) %>%
    group_by(date, SVM_t) %>%
    summarize(N = n()) -> gen_freqs
  
  gen_freqs %>% 
    dcast(SVM_t~date) %>%
    select("CMfall", "CMspring") %>%
    as.matrix() -> Matriz
    Matriz
  
    GTest(Matriz,
        correct="none")
  
  pairwiseNominalIndependence(Matriz,
                              fisher = TRUE,
                              gtest  = TRUE,
                              chisq  = TRUE)
  
 