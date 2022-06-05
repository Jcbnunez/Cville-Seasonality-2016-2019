#### Calcualte overwinter FST across genome
#### 

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
registerDoMC(4)
library(DescTools)
library(scales)


args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])

print(k)

####


####
### Panel a <----------
### FST in Cville over time

# Import a SNP object

SNP_object <- "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST.2.0.poolSNP.Spatial.Temporal.mAF_Miss_Mean_Filt.ECfiltered.Rdata"
load(SNP_object)

filtering_file <- "/project/berglandlab/Dmel_genomic_resources/Filtering_files/snp_dt_25percMissing.Rdata"

inmeta="/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv"
samps <- fread(inmeta)

ingds="/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds"

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

####


########################
### open GDS file
genofile <- seqOpen(ingds, allow.duplicate=T)
#seqClose(genofile)

load(filtering_file)
setkey(snp.dt, id)

use <- apply(snp.dt[,c("VA_ch"), with=F],
             1, any)
snp.dt <- snp.dt[use]
setkey(snp.dt, id)

snp.dt <- snp.dt[N==0 & cm_mb>0 & !is.na(cm_mb) & chr!="X"]
snp.dt %<>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

seqSetFilter(genofile,  variant.id=snp.dt$id)

#g <- snpgdsGetGeno(genofile, snp.id=snp.dt$variant.id)
#snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

## prepare windows
win.bp = 1e6
step.bp = 0.5e6

##setkey(snps.dt , "chr")


wins <- foreach(chr.i=c("2L","2R","3L","3R"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- snp.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)

print(paste(wins$chr[k], wins$start[k], wins$end[k], sep = "-"))


#######

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
snp.dt %>%
  filter(chr == wins$chr[k] & pos > wins$start[k] & pos < wins$end[k] ) %>%
  .$SNP_id ->
  snps_of_win

#for(i in 1:dim(comp_vector_samepops)[1]){

fst.df = foreach(i = 1:dim(comp_vector_samepops)[1], .combine = "rbind" )%dopar%{

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
  
  data.frame(
    samp1 = comp_vector_samepops$samp1[i],
    samp2 = comp_vector_samepops$samp2[i],
    FST = fst.out$FST,
    winSTA = wins$start[k],
    winEND = wins$end[k],
    chr = wins$chr[k]
  )

}## i   


left_join(fst.df, comp_vector_samepops) -> Out_comp_vector_samepops


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
  summarize(mean_f = median(FST))


Filename = paste("obsDat",wins$chr[k], wins$start[k], wins$end[k], "winStats", "Rdata", sep = ".")

save(Out_comp_vector_samepops, 
     file = paste("/scratch/yey2sn/Overwintering_ms/2.Temporal_Spatial_structure/whole_genome_fst_window/",
                  Filename, sep = "") )
