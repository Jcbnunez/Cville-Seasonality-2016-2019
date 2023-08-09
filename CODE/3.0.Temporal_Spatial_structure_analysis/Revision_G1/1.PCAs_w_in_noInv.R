
### Expand the PCA analyses

#vacc
#ijob
#source /gpfs1/home/j/c/jcnunez/.bash_profile
#myspack
#R

### Load packages
library(tidyverse)
library(magrittr)
library(reshape2)
library(FactoMineR)
library(SeqArray)
library(data.table)
library(foreach)
library(scales)

### Load metadata
### Gene annotations
#dmel_annot <- get(load("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/DEST_Cville_SNP_Metadata.Rdata"))

### Load data
load("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/R_EC_filtered_objects/DEST.2.0.poolSNP.Spatial.Temporal.AllDat.ECfiltered.Rdata")
ls()

load("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/R_EC_filtered_objects/DEST_EC_metadata.Rdata")

#### being analysis
#### subset dataset
samps_EFFCOV %>%
  filter(locality %in%     c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")) %>%
  group_by(locality, year) %>%
  summarize(N=n()) %>%
  dcast(locality~year)

samps_EFFCOV %>%
  filter(locality %in%
           c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")) ->
  filtered_samps_for_analysis

####
### Load annotation
annotation <- get(load("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/DEST_Cville_SNP_Metadata.Rdata"))

setDT(annotation)

names(annotation)[1:2] = c("chr","pos")
non.cod = c("intergenic_region","intron_variant","upstream_gene_variant","downstream_gene_variant")
annotation.non.cod = annotation[effect %in% non.cod]

### Annotate inversion
snpdt.obj <- get(load("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/snp_dt.Rdata"))
setDT(snpdt.obj)
snpdt.obj %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
snpdt.obj.NoInv = snpdt.obj[invName == "none"]

### create joint filtering object
dat_filtered_t %>% 
  filter(SNP_id %in% snpdt.obj.NoInv$SNP_id) %>%
  filter(SNP_id %in% annotation.non.cod$SNP_id) ->
  dat_filtered_t.flt

#### 

### quantify

#grep("2L", dat_filtered_t.flt$SNP_id) -> fl2L
#grep("2R", dat_filtered_t.flt$SNP_id) -> fl2R
#grep("3L", dat_filtered_t.flt$SNP_id) -> fl3L
#grep("3R", dat_filtered_t.flt$SNP_id) -> fl3R

which(names(dat_filtered_t.flt) == "SNP_id") -> snpid_colnum
which(names(dat_filtered_t.flt) %in% filtered_samps_for_analysis$sampleId) -> samps.pass


## Some characterizations of AFs and subsequent filtering
MeanAF=c()
MinAF=c()

apply(dat_filtered_t.flt[,-snpid_colnum],
      1, FUN=mean, na.rm=TRUE ) -> MeanAF
data.frame(SNP_id = dat_filtered_t.flt$SNP_id, MeanAF) -> MeanAF

apply(dat_filtered_t.flt,
      1, FUN=min, na.rm=TRUE ) -> MinAF
data.frame(SNP_id = dat_filtered_t.flt$SNP_id, MinAF) -> MinAF

cbind(dat_filtered_t.flt, MeanAF[-1], MinAF[-1]) -> dat_AF_samps_target

##
dat_AF_samps_target %>%
  .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MinAF > 0.001),] ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

###
count_NA = function(x){
  return(sum(is.na(x)))
}
MissDat=c()

pool_cols = grep("SNP_id|MeanAF|MinAF", colnames(dat_AF_samps_target_filtered), invert = T)

apply(dat_AF_samps_target_filtered[,pool_cols],
      1, FUN=count_NA ) -> MissDat

n_pools = length(pool_cols)

data.frame(SNP_id = dat_AF_samps_target_filtered$SNP_id, missing_rate = c(MissDat/n_pools) ) -> MissDat

cbind(dat_AF_samps_target_filtered, MissDat[-1]) -> dat_AF_samps_target_filtered

dat_AF_samps_target_filtered %>%
  filter(missing_rate < 0.01) ->  ### This samples only polymorphic sites
  dat_AF_samps_target_filtered

###

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  summarise(N=n())

dat_AF_samps_target_filtered %>%
  separate(SNP_id, into = c("chr","pos"), sep = "_") %>% 
  group_by(chr) %>%
  summarise(N=n())

#### Create PCAs
load("Fig1.panels.AB.dat.Rdata")
orig.pca = PCA_table

which(names(dat_AF_samps_target_filtered) %in% 
c("SNP_id","MeanAF","MinAF", "missing_rate") ) -> snpid_etc_colnum


oo=
foreach(i = c("2L","2R","3L","3R", "all"),
.combine = "rbind")%do%{

message(i)

if(i != "all"){
grep(i,dat_AF_samps_target_filtered$SNP_id) -> flo
} else if(i == "all"){
flo = dat_AF_samps_target_filtered$SNP_id
}


dat_AF_samps_target_filtered[flo, -snpid_etc_colnum] %>% 
.[, samps.pass] %>% 
t() %>% as.data.frame() %>% 
.[,sample(dim(.)[2], 8000)] %>%
  PCA(graph = F, ncp = 3) ->
  pca.object.flt
data.frame(pca.object.flt$eig) -> eig

pca.object.flt$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleId = rownames(.),
  partitioning = i,
  pc1.var = eig$cumulative.percentage.of.variance[1],
  pc2.var = eig$cumulative.percentage.of.variance[2],
  pc3.var = eig$cumulative.percentage.of.variance[3],
  ) %>% 
  left_join(filtered_samps_for_analysis) -> o

o %>% filter(sampleId %in% orig.pca$sampleId) -> origdat

cor.test(orig.pca$Dim.1, origdat$Dim.1) -> cordat1
cor.test(orig.pca$Dim.2, origdat$Dim.2) -> cordat2
cor.test(orig.pca$Dim.3, origdat$Dim.3) -> cordat3

o %>%
mutate(
corest.1=cordat1$estimate,
corp.1=cordat1$p.value,

corest.2=cordat2$estimate,
corp.2=cordat2$p.value,

corest.3=cordat3$estimate,
corp.3=cordat3$p.value

) -> o.1


return(o.1)

}

### plot
oo %>%
filter(sampleId %in% filtered_samps_for_analysis$sampleId) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    fill = year
  )) +
  geom_point(shape = 21, size = 3.0) +
  facet_grid(~partitioning) +
    scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
  theme_bw() ->
  PCA.grid.dim

ggsave(PCA.grid.dim, file = "PCA.grid.dim.pdf", w = 10, h =2.5)

oo %>%
filter(sampleId %in% filtered_samps_for_analysis$sampleId) %>%
  ggplot(aes(
    x=Dim.2,
    y=Dim.3,
    fill = year
  )) +
  geom_point(shape = 21, size = 3.0) +
  facet_grid(~partitioning) +
    scale_fill_gradientn(
    colors=c("springgreen","cyan","blue","gold","red"),
    values=rescale(c(2011,2013,2015,2016,2018))
  ) +
  theme_bw() ->
  PCA.grid.dim23

ggsave(PCA.grid.dim23, file = "PCA.grid.dim23.pdf", w = 10, h =2.5)

#### bring in data from main PCA
oo %>%
group_by(partitioning) %>%
summarize(
mc1 = mean(corest.1),
mp1 = mean(corp.1),
mc2 = mean(corest.2),
mp2 = mean(corp.2),
mc3 = mean(corest.3),
mp3 = mean(corp.3)
) %>% reshape2::melt(id = "partitioning") %>%
separate(variable, sep = 2, 
into = c("stat", "pc")) %>%
dcast(partitioning+pc~stat, value.var = "value") -> 
cor.values


