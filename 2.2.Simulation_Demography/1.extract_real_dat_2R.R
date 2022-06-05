rm(list = ls())

#load packages
library(adegenet)
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
library(MASS)
library(foreach)
library(doMC)
registerDoMC(2)
library(DescTools)
library(scales)
#install_github('tavareshugo/windowscanr')
#library(windowscanr)

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])

print(k)

euc.dist.2d <- function(coord_vec) {
  #x1, y1, z1, x2, y2, z2
  x1 = coord_vec[1] 
  y1 = coord_vec[2] 
  x2 = coord_vec[3]  
  y2 = coord_vec[4]
  
  sqrt( ((x1 - x2)^2) + ((y1 - y2)^2) )
  
} 

####
### Panel a <----------
### FST in Cville over time

# Import a SNP object

SNP_object <- "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.mAF_Miss_Mean_Filt.ECfiltered.Rdata"
load(SNP_object)

inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
samps <- fread(inmeta)

ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/PoolSNP_noRep_filter/dest.all.PoolSNP.001.50.10Mar2021.ann.noRep.gds"

####

#Generate outfile object
samps %>%
  .[which(.$sampleId %in% rownames(dat_for_Analysis)),] ->
  samps_YwY_EC

samps_YwY_EC$city = gsub("Charlotttesville","Charlottesville", samps_YwY_EC$city)
samps_YwY_EC$city = gsub("Odesa","Odessa", samps_YwY_EC$city )
samps_YwY_EC$city[grep("Yesiloz", samps_YwY_EC$city )] = "Yesiloz"
samps_YwY_EC$city[grep("Kyiv", samps_YwY_EC$city )] = "Kyiv"

print("Modify dates to as.Date format")

samps_YwY_EC %<>% filter(city == "Charlottesville" & year >= 2016)

samps_YwY_EC$collectionDate = as.Date(samps_YwY_EC$collectionDate, 
                                      format='%m/%d/%Y')  

L = dim(samps_YwY_EC)[1]


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
samps_YwY_EC$city %>% unique

##calculate day differences
print("Loop to esrtimate time difference")

for(i in 1:dim(comp_vector)[1]) {
  
  date1=samps_YwY_EC$collectionDate[comp_vector[i,1]]
  date2=samps_YwY_EC$collectionDate[comp_vector[i,2]]
  
  comp_vector$day_diff[i] = abs(as.numeric(date1-date2))
  
  comp_vector$pop1[i] = samps_YwY_EC$city[comp_vector[i,1]]
  comp_vector$pop2[i] = samps_YwY_EC$city[comp_vector[i,2]]
  
  comp_vector$samp1[i] = samps_YwY_EC$sampleId[comp_vector[i,1]]
  comp_vector$samp2[i] = samps_YwY_EC$sampleId[comp_vector[i,2]]
  
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
#seqClose(genofile)
#



### Include DEST sets
samps <- rbind(#samps[set=="DrosRTEC"],
  #samps[set=="DrosEU"],
  samps[set=="CvilleSet"]
)

### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

## choose number of alleles
snps.dt <- snps.dt[nAlleles==2]
snps.dt %<>%
  filter(chr == "2R")
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)
seqSetFilter(genofile, sample.id=samps$sampleId,
             snps.dt[chr%in%c("2R")]$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

## prepare windows
win.bp = 1e6
step.bp = 0.5e6

setkey(snps.dt , "chr")


wins <- foreach(chr.i=c("2R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- snps.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)

print(paste(wins$start[k], wins$end[k], sep = "-"))

#2R	15391154	20276334	2RNS
#wins %<>%
#  filter(start < 1e6 & start > 10e6  )
#dim(wins)


##foreach(i=1:dim(wins)[1], .combine = "rbind")%do%{
##  
##  snps.dt %>%
##    filter(af > 0.05) %>%
##    filter(pos > wins$start[i] & pos < wins$end[i] & chr == wins$chr[i] ) %>%
##    summarize(N = n()) %>%
##    mutate(start = wins$start[i], end = wins$end[i],  chr = wins$chr[i])
##  
##}

#######
snps.dt %<>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))



### select sites

### get allele frequency data
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#dat <- ad/dp
#dim(dat)  
#
### Add metadata
#colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
#rownames(dat) <- seqGetData(genofile, "sample.id")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")


#### filter to windows of interest!
#### filter
snps.dt %>%
  filter( pos > wins$start[k] & pos < wins$end[k] ) %>%
  .$SNP_id ->
  snps_of_win

###################
### Next part
### Calculate FST

#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector_samepops)[1]),
  samp2 = rep(NA, dim(comp_vector_samepops)[1]),
  FST = rep(NA, dim(comp_vector_samepops)[1]),
  winSTA= rep(NA, dim(comp_vector_samepops)[1]),
  winEND=rep(NA, dim(comp_vector_samepops)[1])
)


for(i in 1:dim(comp_vector_samepops)[1]){
  
  print(i/dim(comp_vector_samepops)[1] * 100)
  
  samps_to_compare = c(comp_vector_samepops$samp1[i], comp_vector_samepops$samp2[i])
  
  pool_sizes = c(samps_YwY_EC$nFlies[which(samps_YwY_EC$sampleId == comp_vector_samepops$samp1[i])],
                 samps_YwY_EC$nFlies[which(samps_YwY_EC$sampleId == comp_vector_samepops$samp2[i])])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),snps_of_win]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),snps_of_win]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)
  
  
  fst.out <- computeFST(pool, method = "Anova")
  
  outfile$samp1[i] = comp_vector_samepops$samp1[i]
  outfile$samp2[i] = comp_vector_samepops$samp2[i]
  outfile$FST[i] = fst.out$FST
  outfile$winSTA[i] = wins$start[k]
  outfile$winEND[i] = wins$end[k]
  
}## i   


left_join(comp_vector_samepops, outfile) -> Out_comp_vector_samepops


samps %<>%
  mutate(month_col = month(as.Date(collectionDate, 
                                   format = c("%m/%d/%Y"))))

samps %>%
  dplyr::select(sampleId, month_col, year) -> samp_1_meta
names(samp_1_meta) = c("samp1", "month1", "year1")

samps %>%
  dplyr::select(sampleId, month_col, year) -> samp_2_meta
names(samp_2_meta) =  c("samp2", "month2", "year2")

Out_comp_vector_samepops %<>%
  left_join(samp_1_meta) %>% 
  left_join(samp_2_meta) 

Out_comp_vector_samepops %<>%
  mutate(bin_date = ifelse(.$day_diff <= 200, "1.within", 
                           ifelse(.$day_diff >= 550, "3.Multi-Year", "2.Overwinter" ) ))

Out_comp_vector_samepops %>%
  group_by (bin_date, winSTA, winEND) %>%
  summarize(mean_f = median(FST)) -> FST_values

############## Now do multivariate analysis

dat <- ad/dp
dim(dat)  
#
### Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")


##get metadata
#### load object
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2R.ECfiltered.Rdata")
### ---> filtered_samps_for_analysis

filtered_samps_for_analysis %>%
  filter(MeanEC > 30,
         set == "CvilleSet") %>%
  .$sampleId ->
  samps_for_pca_ld

### generate PCA
dat[samps_for_pca_ld, snps_of_win]  -> dat_raw

dat_raw %>% colMeans() -> means_snps
#names(means_snps[means_snps > 0.05]) %>% length()
names(means_snps[which(means_snps > 0.01)])  -> flt_snps

dat_raw[,flt_snps] %>%
  as.data.frame() %>%
  PCA(graph = F) ->
  PCA_obj_raw

PCA_obj_raw$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(filtered_samps_for_analysis) %>%
  separate(sampleId, remove = F, into = c("pop", "city", "year", "date")) ->
  PCA_dat_tmp

###
###
PCA_obj_raw$eig %>%
  as.data.frame() ->
  variances_explained

### ucl dist
permutations(n = length(PCA_dat_tmp$sampleId), r = 2, repeats.allowed = F, v = PCA_dat_tmp$sampleId) %>%
  as.data.frame()->
  perm_samps

perm_samps %>%
  dplyr::select(sampleId = V1) %>%
  left_join(PCA_dat_tmp[,c( "sampleId","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s1 = sampleId,
                Dim.1.s1 = Dim.1,
                Dim.2.s1 = Dim.2,
                Dim.3.s1 =Dim.3
  ) -> S1_dat

perm_samps %>%
  dplyr::select(sampleId = V2) %>%
  left_join(PCA_dat_tmp[,c( "sampleId","Dim.1", "Dim.2", "Dim.3")]) %>%
  dplyr::select(sampleId.s2 = sampleId,
                Dim.1.s2 = Dim.1,
                Dim.2.s2 = Dim.2,
                Dim.3.s2 =Dim.3
  ) -> S2_dat

cbind(S1_dat, S2_dat) -> euc.in.dat

euc.in.dat %>%
  mutate(Dim.1.s1.mean = mean(Dim.1.s1),
         Dim.2.s1.mean = mean(Dim.2.s1),
         Dim.1.s2.mean = mean(Dim.1.s2),
         Dim.2.s2.mean = mean(Dim.2.s2),
         Dim.1.s1.min = min(Dim.1.s1),
         Dim.2.s1.min = min(Dim.2.s1),
         Dim.1.s2.min = min(Dim.1.s2),
         Dim.2.s2.min = min(Dim.2.s2),
         Dim.1.s1.max = max(Dim.1.s1),
         Dim.2.s1.max = max(Dim.2.s1),
         Dim.1.s2.max = max(Dim.1.s2),
         Dim.2.s2.max = max(Dim.2.s2)
  ) %>%
  mutate(
    Dim.1.s1.STD = (Dim.1.s1-Dim.1.s1.mean)/(Dim.1.s1.max-Dim.1.s1.min), 
    Dim.2.s1.STD = (Dim.2.s1-Dim.2.s1.mean)/(Dim.2.s1.max-Dim.2.s1.min), 
    Dim.1.s2.STD = (Dim.1.s2-Dim.1.s2.mean)/(Dim.1.s2.max-Dim.1.s2.min), 
    Dim.2.s2.STD = (Dim.2.s2-Dim.2.s2.mean)/(Dim.2.s2.max-Dim.2.s2.min)
  ) -> euc.in.dat.std


eu_dist.2d.D.pca = foreach(i = 1:nrow(euc.in.dat.std), .combine = c )%do%{
  
  in_coords = unlist(c(euc.in.dat.std[i, c("Dim.1.s1.STD",  "Dim.2.s1.STD")],
                       euc.in.dat.std[i, c("Dim.1.s2.STD",  "Dim.2.s2.STD")]))
  
  euc.dist.2d(as.vector(in_coords))
  
}

cbind(euc.in.dat.std, eu_dist.2d.D.pca) %>%
  as.data.frame %>% 
  separate(sampleId.s1, into = c("pop1", "city1", "year1", "date1"), sep = "_", remove = F) %>%
  separate(sampleId.s2, into = c("pop2", "city2", "year2", "date2"), sep = "_", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) ->
  Euclidean.analysis.results.2d.pca


Euclidean.analysis.results.2d.pca %>%
  group_by(as.factor(year_diff)) %>%
  summarize(Median.standarized = quantile(eu_dist.2d.D.pca, 0.5),
            IQR05.standarized = quantile(eu_dist.2d.D.pca, 0.05),
            IQR95.standarized = quantile(eu_dist.2d.D.pca, 0.95)
  ) -> EUC_PCA_standarized

### Do correlation PCA
PCA_dat_tmp %>% 
  dplyr::select(Dim.1, Dim.2, Dim.3, year, MeanEC, `In(2L)t`, `In(2R)Ns`) %>%
  melt(id = c("year", "MeanEC", "In(2L)t", "In(2R)Ns" ) ) %>%
  mutate(year = as.numeric(year)) %>%
  group_by(variable) %>%
  summarize(corr.time.est = cor.test(value, year)$est,
            corr.time.p = cor.test(value, year)$p.value,
            corr.EC.est = cor.test(value, MeanEC)$est,
            corr.EC.p = cor.test(value, MeanEC)$p.value,
            corr.2lt.est = cor.test(value, `In(2L)t`)$est,
            corr.2lt.p = cor.test(value, `In(2L)t`)$p.value,
            corr.2rns.est = cor.test(value, `In(2R)Ns`)$est,
            corr.2rns.p = cor.test(value, `In(2R)Ns`)$p.value,
  ) -> correlations_tmp

### Do correlation DACP
dat_raw[,flt_snps] %>%
  as.data.frame() %>%
  dapc(. , grp = as.factor(filter(filtered_samps_for_analysis, MeanEC >= 30, locality == "VA_ch", year >= 2016 )$year), 
       n.pca=9, n.da=5 ) ->
  dapc_first_pass

dat_raw[,flt_snps] %>%
  as.data.frame() %>%
  dapc(. , grp = as.factor(filter(filtered_samps_for_analysis, MeanEC >= 30, locality == "VA_ch", year >= 2016 )$year),
       n.pca=optim.a.score(dapc_first_pass)$best, n.da=5 ) ->
  dapc_optim

dapc_optim$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(filtered_samps_for_analysis) ->
  dapc_optim_dat

## Correlations models
dapc_optim_dat %>%
  dplyr::select(LD1,LD2, year, MeanEC) %>%
  melt(id = c("year","MeanEC") ) ->
  dapc_table_melt

dapc_table_melt  %>% 
  group_by(variable) %>%
  summarize(corr.time.est = cor.test(value, year)$est,
            corr.time.p = cor.test(value, year)$p.value,
            corr.EC.est = cor.test(value, MeanEC)$est,
            corr.EC.p = cor.test(value, MeanEC)$p.value,
  ) -> ld_correlations


## Euclidean distance
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


euc.in.dat.dapc %>%
  mutate(LD1.s1.mean = mean(LD1.s1),
         LD2.s1.mean = mean(LD2.s1),
         LD1.s2.mean = mean(LD1.s2),
         LD2.s2.mean = mean(LD2.s2),
         LD1.s1.min = min(LD1.s1),
         LD2.s1.min = min(LD2.s1),
         LD1.s2.min = min(LD1.s2),
         LD2.s2.min = min(LD2.s2),
         LD1.s1.max = max(LD1.s1),
         LD2.s1.max = max(LD2.s1),
         LD1.s2.max = max(LD1.s2),
         LD2.s2.max = max(LD2.s2)
  ) %>%
  mutate(
    LD1.s1.STD = (LD1.s1-LD1.s1.mean)/(LD1.s1.max-LD1.s1.min), 
    LD2.s1.STD = (LD2.s1-LD2.s1.mean)/(LD2.s1.max-LD2.s1.min), 
    LD1.s2.STD = (LD1.s2-LD1.s2.mean)/(LD1.s2.max-LD1.s2.min), 
    LD2.s2.STD = (LD2.s2-LD2.s2.mean)/(LD2.s2.max-LD2.s2.min)
  ) -> euc.in.dat.dapc.std


eu_dist.2d.D = foreach(i = 1:nrow(euc.in.dat.dapc.std), .combine = c )%do%{
  
  in_coords = unlist(c(euc.in.dat.dapc.std[i, c("LD1.s1",  "LD2.s1")],
                       euc.in.dat.dapc.std[i, c("LD1.s2",  "LD2.s2")]))
  
  euc.dist.2d(as.vector(in_coords))
  
}

cbind(euc.in.dat.dapc.std, eu_dist.2d.D) %>%
  as.data.frame %>% 
  separate(sampleId.s1, into = c("pop1", "city1", "year1", "date1"), sep = "_", remove = F) %>%
  separate(sampleId.s2, into = c("pop2", "city2", "year2", "date2"), sep = "_", remove = F) %>% 
  mutate(year_diff = abs(as.numeric(year1) - as.numeric(year2))) ->
  Euclidean.analysis.results.dapc

Euclidean.analysis.results.dapc %>%
  group_by(as.factor(year_diff)) %>%
  summarize(Median.standarized = quantile(eu_dist.2d.D, 0.5),
            IQR05.standarized = quantile(eu_dist.2d.D, 0.05),
            IQR95.standarized = quantile(eu_dist.2d.D, 0.95))  -> EUC_DAPC_standarized

### Variance of Allele frequency
dat_raw[,flt_snps] %>%
  as.data.frame() %>%
  apply(., 2, sd) %>% quantile(c(0.05, 0.25, 0.5, 0.75, 0.95)) -> variance_in_af



##### outputs!
##PCA
correlations_tmp %>%
  mutate(corr.time.est = abs(corr.time.est)) %>%
  dplyr::select(corr.time.est, Dim=variable) %>%
  melt() %>%
  mutate(Obs_var = paste(paste("Dim=", Dim, sep = ""),variable, sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt1

variances_explained[1:3,] %>%
  dplyr::select( "percentage of variance"  ) %>%
  mutate(Dim = rownames(.)) %>%
  melt  %>%
  mutate(Obs_var = paste(paste("DimPCA=", Dim, sep = ""),variable, sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt2
pt2$Obs_var = gsub( " ", "_" , pt2$Obs_var)

EUC_PCA_standarized %>%
  melt %>%
  mutate(Obs_var = paste(paste("yearDiff=", `as.factor(year_diff)`, sep = ""), 
                         paste(variable, "PCAEucDist", sep  = ".") , sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt3

## DAPC
ld_correlations %>%
  mutate(corr.time.est = abs(corr.time.est)) %>%
  dplyr::select(corr.time.est, Dim=variable) %>%
  melt() %>%
  mutate(Obs_var = paste(paste("DimDAPC=", Dim, sep = ""),variable, sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt4


EUC_DAPC_standarized %>%
  melt %>%
  mutate(Obs_var = paste(paste("yearDiff=", `as.factor(year_diff)`, sep = ""), 
                         paste(variable, "DAPCEucDist", sep  = ".") , sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt5

## Variance in allele frequency
variance_in_af %>%
  data.frame(AF_var=.) %>%
  mutate(Var_AF_quantile = rownames(.)) %>%
  melt %>%
  mutate(Obs_var = paste(paste("Global=", 
                               Var_AF_quantile, sep = ""), variable, sep = "_")) %>%
  dplyr::select(Obs_var, value) -> pt6


### FST values
FST_values %>%
  .[,c("bin_date", "mean_f")] %>%
  mutate(Obs_var = paste(paste("yearDiff=", 
                               bin_date, sep = ""), "FST", sep = "_")) %>%
  melt %>%
  dplyr::select(Obs_var, value) -> pt7

#### Merge datasets  
rbind(pt1, pt2, pt3, pt4, pt5, pt6, pt7) %>%
  as.data.frame() %>%
  mutate(win.start = wins$start[k],
         win.end = wins$end[k]) ->
  data.out

Filename = paste("obsDat", wins$start[k], wins$end[k], "winStats", "Rdata", sep = ".")

save(data.out, file = Filename)

