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
library(adegenet)
library(scales)
library(lme4)
library(foreach)

######
## define Euclidean distances
euc.dist.3d <- function(coord_vec) {
  #x1, y1, z1, x2, y2, z2
  x1 = coord_vec[1] 
  y1 = coord_vec[2] 
  z1 = coord_vec[3] 
  x2 = coord_vec[4]  
  y2 = coord_vec[5]
  z2 = coord_vec[6]
  
  sqrt( ((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))
  
} 

## define Euclidean distances
euc.dist.2d <- function(coord_vec) {
  #x1, y1, z1, x2, y2, z2
  x1 = coord_vec[1] 
  y1 = coord_vec[2] 
  x2 = coord_vec[3]  
  y2 = coord_vec[4]

  sqrt( ((x1 - x2)^2) + ((y1 - y2)^2) )
  
} 


#####
#####
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
           MeanEC > 30,
           year >= 2016 ) %>% 
    .$sampleId -> select_samples
  
  o %>%
    as.data.frame() %>% 
    filter(rownames(.) %in%  select_samples) ->
    snps.tmp
  
  pca_table_list[[i]] = snps.tmp
  
}

chrs_o = do.call(cbind, pca_table_list)

### Calculate variances of allele frequency
### 
chrs_o %>% colnames -> SNP_ids
chrs_o %>% rownames -> samp_ids

chrs_o %>%
  t %>%
  as.data.frame() ->
  chrs_o_t

##calculate the variance per site
apply(chrs_o_t, 1, var) -> VAR_per_site

VAR_per_site %>% mean
#0.004244469

### build PCA
### 
chrs_o %>%
  as.data.frame() %>% 
  filter(rownames(.) %in%  select_samples) ->
  chrs_o_flt
  
rownames(chrs_o_flt) %>%
  data.frame(sampleId = .) %>%
  separate(sampleId, into = c("pop", "city", "year", "season"), sep = "_"  ) ->
  samps_metadat

####
####
####
#### Do PCA

chrs_o_flt %>%
  PCA(scale.unit = F, graph = F, ncp = 20) ->
  PCA_object

PCA_object$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.),
         chr = lead_snp_guide$CHR) %>%  
  left_join(., filtered_samps_for_analysis ) -> PCA_table

## PC 12
PCA_table %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill=year,
  )) +
  geom_point( size = 4, shape = 21) +
  scale_fill_gradientn(
    colors=c("blue","gold","red"),
    values=rescale(c(2016,2017,2018))
  ) -> pca_plot12

ggsave(pca_plot12, file = "pca_plot12.pdf", w= 5, h =4)

### PC 23
PCA_table %>%
  ggplot(aes(
    x=Dim.2,
    y=Dim.3,
    fill=year,
  )) +
  geom_point( size = 4, shape = 21) +
  scale_fill_gradientn(
    colors=c("blue","gold","red"),
    values=rescale(c(2016,2017,2018))
  ) -> pca_plot23

ggsave(pca_plot23, file = "pca_plot23.pdf", w= 5, h =4)

### Some analysis ### Some analysis

## Linear models
PCA_table %>%
  dplyr::select(Dim.1,Dim.2,Dim.3, year, MeanEC) %>%
  melt(id = c("year","MeanEC") ) ->
  PCA_table_melt
lmList(data =PCA_table_melt, value ~ year + MeanEC | variable ) %>% summary()

## Correlations models
PCA_table_melt  %>% 
  group_by(variable) %>%
  summarize(corr.time.est = cor.test(value, year)$est,
            corr.time.p = cor.test(value, year)$p.value,
            corr.EC.est = cor.test(value, MeanEC)$est,
            corr.EC.p = cor.test(value, MeanEC)$p.value,
            )

## Euclidean distances
## define comparions
## library(gtools) ## for the permutation commands

permutations(n = length(PCA_table$sampleId), r = 2, repeats.allowed = F, v = PCA_table$sampleId) %>%
  as.data.frame()->
  perm_samps

perm_samps %>%
  dplyr::select(sampleId = V1) %>%
  left_join(PCA_table[,c( "sampleId","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s1 = sampleId,
                Dim.1.s1 = Dim.1,
                Dim.2.s1 = Dim.2,
                Dim.3.s1 =Dim.3
                ) -> S1_dat
  
perm_samps %>%
  dplyr::select(sampleId = V2) %>%
  left_join(PCA_table[,c( "sampleId","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s2 = sampleId,
                Dim.1.s2 = Dim.1,
                Dim.2.s2 = Dim.2,
                Dim.3.s2 =Dim.3
  ) -> S2_dat

cbind(S1_dat, S2_dat) -> euc.in.dat


eu_dist.3d.D = foreach(i = 1:nrow(euc.in.dat), .combine = c )%do%{
  
  in_coords = unlist(c(euc.in.dat[i, c("Dim.1.s1",  "Dim.2.s1", "Dim.3.s1")],
                       euc.in.dat[i, c("Dim.1.s2",  "Dim.2.s2", "Dim.3.s2")]))
  
  euc.dist.3d(as.vector(in_coords))

}

eu_dist.2d.D.pca = foreach(i = 1:nrow(euc.in.dat), .combine = c )%do%{
  
  in_coords = unlist(c(euc.in.dat[i, c("Dim.1.s1",  "Dim.2.s1")],
                       euc.in.dat[i, c("Dim.1.s2",  "Dim.2.s2")]))
  
  euc.dist.2d(as.vector(in_coords))
  
}

  
cbind(euc.in.dat, eu_dist.3d.D) %>%
  as.data.frame %>% 
  separate(sampleId.s1, into = c("pop1", "city1", "year1", "date1"), sep = "_", remove = F) %>%
  separate(sampleId.s2, into = c("pop2", "city2", "year2", "date2"), sep = "_", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) ->
  Euclidean.analysis.results

Euclidean.analysis.results %>%
  ggplot(aes(
    x=as.factor(year_diff),
    y=eu_dist.3d.D
  )) + 
  geom_boxplot() ->
  euc_box

ggsave(euc_box, file ="euc_box.pdf", w = 4, h = 4)

lm(eu_dist.3d.D~year_diff, data = Euclidean.analysis.results) %>% summary()

Euclidean.analysis.results %>%
  group_by(as.factor(year_diff)) %>%
  summarize(Median = quantile(scale(eu_dist.3d.D), 0.5),
            IQR05 = quantile(scale(eu_dist.3d.D), 0.05),
            IQR95 = quantile(scale(eu_dist.3d.D), 0.95)
            )


cbind(euc.in.dat, eu_dist.2d.D.pca) %>%
  as.data.frame %>% 
  separate(sampleId.s1, into = c("pop1", "city1", "year1", "date1"), sep = "_", remove = F) %>%
  separate(sampleId.s2, into = c("pop2", "city2", "year2", "date2"), sep = "_", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) ->
  Euclidean.analysis.results.2d.pca

Euclidean.analysis.results.2d.pca %>%
  ggplot(aes(
    x=as.factor(year_diff),
    y=eu_dist.2d.D.pca
  )) + 
  geom_boxplot() ->
  euc_box.2d.pca
ggsave(euc_box.2d.pca, file ="euc_box.2d.pca.pdf", w = 4, h = 4)

lm(eu_dist.2d.D.pca~year_diff, data = Euclidean.analysis.results.2d.pca) %>% summary()

Euclidean.analysis.results %>%
  group_by(as.factor(year_diff)) %>%
  summarize(Median = quantile(eu_dist.3d.D, 0.5),
            IQR05 = quantile(eu_dist.3d.D, 0.05),
            IQR95 = quantile(eu_dist.3d.D, 0.95)
  )



####
####
####
#### Do DAPC
chrs_o_flt %>%
  dapc(. , grp = as.factor(samps_metadat$year), n.pca=9, n.da=5 ) ->
  dapc_first_pass

chrs_o_flt %>%
  dapc(. , grp = as.factor(samps_metadat$year), n.pca=optim.a.score(dapc_first_pass)$best, n.da=5 ) ->
  dapc_optim

pdf("alpha.score.pdf")
optim.a.score(dapc_first_pass)
dev.off()

dapc_optim$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(filtered_samps_for_analysis) ->
  dapc_optim_dat

###
## dapc 12
dapc_optim_dat %>%
  ggplot(aes(
    x=LD1,
    y=LD2,
    fill=year,
  )) +
  geom_point( size = 4, shape = 21) +
  scale_fill_gradientn(
    colors=c("blue","gold","red"),
    values=rescale(c(2016,2017,2018))
  ) -> dapc_plot12

ggsave(dapc_plot12, file = "dapc_plot12.pdf", w= 5, h =4)

## Linear models
dapc_optim_dat %>%
  dplyr::select(LD1,LD2, year, MeanEC) %>%
  melt(id = c("year","MeanEC") ) ->
  dapc_table_melt

lmList(data =dapc_table_melt, value ~ year + MeanEC | variable ) %>% summary()
lmList(data =dapc_table_melt, value ~ year + MeanEC | variable ) %>% summary() %>% .$r.squared

## Correlations models
dapc_table_melt  %>% 
  group_by(variable) %>%
  summarize(corr.time.est = cor.test(value, year)$est,
            corr.time.p = cor.test(value, year)$p.value,
            corr.EC.est = cor.test(value, MeanEC)$est,
            corr.EC.p = cor.test(value, MeanEC)$p.value,
  )


##Euclidean distance

permutations(n = length(dapc_optim_dat$sampleId), r = 2, repeats.allowed = F, v = dapc_optim_dat$sampleId) %>%
  as.data.frame()->
  perm_samps.dapc

perm_samps.dapc %>%
  dplyr::select(sampleId = V1) %>%
  left_join(dapc_optim_dat[,c( "sampleId","LD1", "LD2")]) %>%
  dplyr::select(sampleId.s1 = sampleId,
                LD1.s1 = LD1,
                LD2.s1 = LD2
  ) -> S1_dat.dapc

perm_samps.dapc %>%
  dplyr::select(sampleId = V2) %>%
  left_join(dapc_optim_dat[,c( "sampleId","LD1", "LD2")]) %>%
  dplyr::select(sampleId.s2 = sampleId,
                LD1.s2 = LD1,
                LD2.s2 = LD2
  ) -> S2_dat.dapc

cbind(S1_dat.dapc, S2_dat.dapc) -> euc.in.dat.dapc


eu_dist.2d.D = foreach(i = 1:nrow(euc.in.dat.dapc), .combine = c )%do%{
  
  in_coords = unlist(c(euc.in.dat.dapc[i, c("LD1.s1",  "LD2.s1")],
                       euc.in.dat.dapc[i, c("LD1.s2",  "LD2.s2")]))
  
  euc.dist.2d(as.vector(in_coords))
  
}


cbind(euc.in.dat.dapc, eu_dist.2d.D) %>%
  as.data.frame %>% 
  separate(sampleId.s1, into = c("pop1", "city1", "year1", "date1"), sep = "_", remove = F) %>%
  separate(sampleId.s2, into = c("pop2", "city2", "year2", "date2"), sep = "_", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) ->
  Euclidean.analysis.results.dapc

Euclidean.analysis.results.dapc %>%
  ggplot(aes(
    x=as.factor(year_diff),
    y=eu_dist.2d.D
  )) + 
  geom_boxplot() ->
  euc_box.dapc

ggsave(euc_box.dapc, file ="euc_box.dapc.pdf", w = 4, h = 4)

lm(eu_dist.2d.D~year_diff, data = Euclidean.analysis.results.dapc) %>% summary()

Euclidean.analysis.results.dapc %>%
  group_by(as.factor(year_diff)) %>%
  summarize(Median = quantile(eu_dist.2d.D, 0.5),
            IQR05 = quantile(eu_dist.2d.D, 0.05),
            IQR95 = quantile(eu_dist.2d.D, 0.95)
  )
