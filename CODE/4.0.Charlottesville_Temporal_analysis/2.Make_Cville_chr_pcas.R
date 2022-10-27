### make PCA from Cville chromosomes 
### 

rm(list = ls())
# Load packages

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
#install_github('tavareshugo/windowscanr')
#library(windowscanr)

######
###### Import R objects

objects <- c(
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2R.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3L.ECfiltered.Rdata",
"/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3R.ECfiltered.Rdata"
)

pca_table_list = list()
for(i in 1:length(objects)){
  
load(objects[i])
o %>% colnames() %>% .[1] -> lead_snp

lead_snp %>% data.frame(headsnp = .) %>% 
  separate(headsnp , into = c("CHR","POS")) -> lead_snp_guide

print(lead_snp_guide$CHR)

filtered_samps_for_analysis %>%
  filter(city == "Charlottesville",
         MeanEC > 30) %>% 
  .$sampleId -> select_samples

o %>%
  as.data.frame() %>% 
  filter(rownames(.) %in%  select_samples) %>% 
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  PCA_object

PCA_object$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.),
         chr = lead_snp_guide$CHR) %>%  
  left_join(., filtered_samps_for_analysis ) -> PCA_table

pca_table_list[[i]] = PCA_table
}

pca_table_df = do.call(rbind, pca_table_list)

### add some extra time metadata

pca_table_df %<>%
mutate(month_col = month(as.Date(collectionDate, 
                                 format = c("%m/%d/%Y"))))
  

save(pca_table_df, file = "make.Fig2.panelA.Rdata")
#### general plot
#### Make the arrows


  ggplot() + 
  geom_point(shape =21, 
             size = 3,
             alpha = 0.9,
             data = pca_table_df,
             aes(
               x=Dim.1,
               y=Dim.2,
               fill = year)
             ) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "gold",
                       midpoint = 2016) +
  theme_bw() +
    ylab("PC2") +
    xlab("PC1") +
  facet_wrap(~chr) ->
  cville_chr_pca

ggsave(cville_chr_pca,
       file = "cville_chr_pca.pdf",
       width = 5,
       height = 4)

##########
### plot specific panels
#####

## inv 2LT
pca_table_df %>%
  ggplot(aes(
    x=`In(2L)t`,
    y=Dim.1,
    color = chr
  )) + 
  geom_smooth(method = "lm", se = F) +
  geom_point(size = 2)  ->
  inv_vs_pc1

ggsave(inv_vs_pc1,
       file = "inv_vs_pc1.pdf",
       width = 4,
       height = 4)

## year signal
pca_table_df %>%
  ggplot(aes(
    x=year,
    y=Dim.2,
    color = chr
  )) + 
  geom_smooth(method = "lm", se = F) +
  geom_point(size = 2)  ->
  year_vs_pc2

ggsave(year_vs_pc2,
       file = "year_vs_pc2.pdf",
       width = 4,
       height = 4)

## year signal
pca_table_df %>%
  ggplot(aes(
    x=as.factor(month_col),
    y=abs(Dim.1),
    fill = chr
  )) + 
  geom_boxplot()  ->
  month_vs_pc1

ggsave(
inv_vs_pc1/
year_vs_pc2/
month_vs_pc1,
       file = "pannel_figures.pdf",
       width = 4,
       height = 9)

########################
########################
##### RUN REGRESSIONS

###Correlations with 2lt? -- begin
PCA_object$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.),
         Run = chr) %>% 
  left_join(., filtered_samps_for_analysis ) -> PCA_table

dims_captured=dim(PCA_object$ind$coord)[2]

PCA_table$year = as.numeric(PCA_table$year)

### Save 
save(PCA_table, file = paste(chr, "pca_coords", "Rdata", sep = "."))

#### Run linear regressions
#### 
guide=paste("Dim", 1:dims_captured, sep = ".")
Dim_lm_list = list()

PCA_table %>%
  .[,c(grep("Dim", colnames(.)),
       grep("year", colnames(.)))] %>%
  melt(id = "year") ->
  collection_year_df

PCA_model_fit = collection_year_df
type="year"

for(i in 1:dims_captured){
  print(i)
  
  lm(data = PCA_model_fit[which(PCA_model_fit$variable == guide[i] ),],
     value ~ year ) -> tmp
  
  summary(tmp) %>%
    .$coefficients %>%
    as.data.frame() %>%
    mutate(Dim = guide[i], 
           terms = rownames(.),
           run_type= chr,
           test_type = type
    ) -> tmp_s
  
  Dim_lm_list[[i]] = tmp_s[grep("(Intercept)", rownames(tmp_s), invert = T),]
  
}

Dim_lm_time = do.call(rbind, Dim_lm_list)
Dim_lm_time %<>%
  mutate(p.val.fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
         p.val.bonfe = p.adjust(`Pr(>|t|)`, method = "bonferroni"))

############# Inversions
############# 

inversions =c(
  "In(2L)t",
  "In(2R)Ns",
  "In(3L)P",
  "In(3R)C",
  "In(3R)Mo",
  "In(3R)Payne")

inversions = gsub("[\\(\\)]", "_", inversions)
names(PCA_table) = gsub("[\\(\\)]", "_", names(PCA_table))


PCA_table %>%
  .[,c(grep("Dim", colnames(.)),which(names(.) %in% inversions))] %>%
  melt(id = inversions) ->
  in_df


guide=paste("Dim", 1:dims_captured, sep = ".")
Dim_lm_list = list()
PCA_model_fit = in_df
type="inv"
for(i in 1:dims_captured){
  print(i)
  
  
  varlist <- inversions
  
  models <- lapply(varlist, function(x) {
    lm(substitute(value ~ i, list(i = as.name(x))), 
       data = PCA_model_fit[which(PCA_model_fit$variable == guide[i] ),])
  })
  
  lapply(models, summary) %>%
    lapply(., coef) %>%
    do.call(rbind, .) %>%
    .[grep("(Intercept)", rownames(.), invert = T),] %>%
    as.data.frame() %>%
    mutate(Dim = guide[i], 
           terms = rownames(.),
           run_type= chr,
           test_type = type) -> tmp_s
  
  Dim_lm_list[[i]] = tmp_s[]
  
}

Dim_lm_inv = do.call(rbind, Dim_lm_list)
Dim_lm_inv %<>%
  mutate(p.val.fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
         p.val.bonfe = p.adjust(`Pr(>|t|)`, method = "bonferroni"))

############# Effective Coverage
############# 

PCA_table %>%
  .[,c(grep("Dim", colnames(.)),
       grep("MeanEC", colnames(.)))] %>%
  melt(id = "MeanEC") ->
  MeanEC_year_df

guide=paste("Dim", 1:dims_captured, sep = ".")
Dim_lm_list = list()
PCA_model_fit = MeanEC_year_df
type="EC"
for(i in 1:dims_captured){
  print(i)
  
  lm(data = PCA_model_fit[which(PCA_model_fit$variable == guide[i] ),],
     value ~ MeanEC ) -> tmp
  
  summary(tmp) %>%
    .$coefficients %>%
    as.data.frame() %>%
    mutate(Dim = guide[i], 
           terms = rownames(.),
           run_type= chr,
           test_type = type
    ) -> tmp_s
  
  Dim_lm_list[[i]] = tmp_s[grep("(Intercept)", rownames(tmp_s), invert = T),]
  
}

Dim_lm_EC = do.call(rbind, Dim_lm_list)
Dim_lm_EC %<>%
  mutate(p.val.fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
         p.val.bonfe = p.adjust(`Pr(>|t|)`, method = "bonferroni"))

#### save output
#### 

rbind(
  Dim_lm_time,
  Dim_lm_inv,
  Dim_lm_EC
) -> merged_lm_output

write.table(merged_lm_output, 
            file = paste(chr, "lm_fits", "txt" , sep = ".") ,
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = TRUE)

