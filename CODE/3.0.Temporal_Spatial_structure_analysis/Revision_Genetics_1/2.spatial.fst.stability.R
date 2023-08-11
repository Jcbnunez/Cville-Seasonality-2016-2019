rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
kor=as.numeric(args[1])
### Shift is a parameter to add to k, to compensate for the array limit of VACC
shift=as.numeric(args[2])

print(kor)
print(shift)

k=kor+shift
print(k)

library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(SeqArray)
library(lubridate)
library(tidyverse)
library(gtools)
library(poolfstat)
library(geosphere)


#### load snp metadata

#### import genomic data
#cp /netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/dest.all.PoolSNP.001.50.10Mar2021.ann.gds ./

#cp /netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/R_EC_filtered_objects/DEST_EC_metadata.Rdata ./

#cp /netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/DEST_Cville_SNP_Metadata.Rdata ./

#cp /netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/snp_dt.Rdata ./

genofile <- seqOpen("../dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

#### load metadata
load("../DEST_EC_metadata.Rdata")
#samps_EFFCOV %>%
#  filter(locality %in%     c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")) %>%
#  group_by(locality, year) %>%
#  summarize(N=n()) %>%
#  dcast(locality~year)
# --> generate the samps file

temporal_samples = c("DE_Bro","DE_Mun","FI_Aka","PA_li","TR_Yes","UA_Ode", "UA_od","VA_ch","WI_cp")

samps_EFFCOV %>% head
samps = samps_EFFCOV %>% filter(MeanEC >= 28)
###
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")

snps.dt <- snps.dt[nAlleles==2][missing < 0.1][chr %in% c("2L","2R","3L","3R")]

snps.dt %<>% mutate(SNP_id = paste(chr, pos, sep = "_")) 
snps.dt %>% dim

### Load annotation
annotation <- get(load("DEST_Cville_SNP_Metadata.Rdata"))

setDT(annotation)

names(annotation)[1:2] = c("chr","pos")
non.cod = c("intergenic_region","intron_variant","upstream_gene_variant","downstream_gene_variant")
annotation.non.cod = annotation[effect %in% non.cod]


### Annotate inversion
snpdt.obj <- get(load("snp_dt.Rdata"))
setDT(snpdt.obj)
snpdt.obj %<>% mutate(SNP_id = paste(chr, pos, sep = "_"))
snpdt.obj.NoInv = snpdt.obj[invName == "none"]

#####
###
snps.dt %>% 
  filter(SNP_id %in% snpdt.obj.NoInv$SNP_id) %>%
  filter(SNP_id %in% annotation.non.cod$SNP_id) ->
  snps.dt.FLT

###
snps.dt %>% 
  filter(SNP_id %in% snpdt.obj.NoInv$SNP_id) %>%
  filter(SNP_id %in% annotation.non.cod$SNP_id) ->
  snps.dt.FLT

### Exclude samples that are other sp. or far-away
samps %>%
group_by(year) %>%
summarize(N = n())


samps %>%
group_by(locality, year) %>%
slice_head %>%
group_by(year) %>%
summarize(N = n())

samps %>%
filter(year %in% 2014:2016) ->
samps.2014.2016

samps.2014.2016 %>%
group_by(locality, year) %>%
slice_head %>%
group_by(continent) %>%
summarize(N = n())

### --- extract data ---- 
seqSetFilter(genofile, 
sample.id=samps.2014.2016$sampleId,
variant.id=snps.dt.FLT$variant.id)

ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

sampleids <- seqGetData(genofile, "sample.id")

colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")



### create the data object
#dat = ad/dp

#### Create matrix of comparisons

year.fst.mat =
foreach(y = 2014:2016, .combine = "rbind")%do%{

tmp.samps = samps.2014.2016 %>%
filter(year == y) %>%
group_by(locality) %>%
slice_head

message(y)

tmp.samps$collectionDate = as.Date(tmp.samps$collectionDate, 
                                      format='%m/%d/%Y')  

L = dim(tmp.samps)[1]

comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

comp.vector =
foreach(i = 1:dim(comp_vector)[1],
.combine = "rbind")%do%{

  message(i)

  samp1 = tmp.samps$sampleId[comp_vector[i,1]]
  samp2 = tmp.samps$sampleId[comp_vector[i,2]]
 
  lat.1 = tmp.samps$lat[comp_vector[i,1]]
  long.1 = tmp.samps$long[comp_vector[i,1]]
  
  lat.2 = tmp.samps$lat[comp_vector[i,2]]
  long.2 = tmp.samps$long[comp_vector[i,2]]
  
  pool1 = tmp.samps$nFlies[comp_vector[i,1]]
  pool2 = tmp.samps$nFlies[comp_vector[i,2]]
 
  cont1 = tmp.samps$continent[comp_vector[i,1]]
  cont2 = tmp.samps$continent[comp_vector[i,2]]

 
  data.frame(
  samp1,
  samp2,
  cont1,
  cont2,
  lat.1 ,
  long.1,
  lat.2,
  long.2,
  pool1,
  pool2,
  year=y
  ) %>%
   mutate(hav_d = distHaversine(
        matrix(c(long.1, long.2, lat.1, lat.2),
        nrow = 2))/1000)
  
} ## i

return(comp.vector)

}

year.fst.mat %>%
filter(cont1==cont2) ->
year.fst.mat

year.fst.mat %>%
group_by(cont1, cont2 ) %>%
summarize(N = n())

#### --> launch analysis
#year.fst.mat

#o =
#foreach(
#k = 1:dim(year.fst.mat)[1] #,
# .combine = "rbind",
# .errorhandling = "remove")%do%{

#print(k/dim(year.fst.mat)[1] * 100)

#k = 1
##### ===>>> K is declared above

  samps_to_compare = c(year.fst.mat$samp1[k], year.fst.mat$samp2[k])
  pool_sizes = c(year.fst.mat$pool1[k], year.fst.mat$pool2[k])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  

  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2)

  fst.out <- computeFST(pool, method = "Anova")
  
  year.fst.mat[k,] %>%
  mutate(FST = fst.out$FST) ->
  o.tmp
 
# return(o.tmp)

#}
root = "/gpfs2/scratch/jcnunez/genetics_resub/2.spatial_fst/spatial.fst/"

save(o.tmp,
file = paste(root, "spacfst." ,k, ".Rdata", sep = ""))
