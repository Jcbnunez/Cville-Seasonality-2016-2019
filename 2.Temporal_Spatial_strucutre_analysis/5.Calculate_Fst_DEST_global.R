rm(list = ls())

#load packages
library(tidyverse)
library(magrittr)
library(forcats)
library(FactoMineR)
library(slider)
library(gtools)
library(poolfstat)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gmodels)
#install_github('tavareshugo/windowscanr')
library(windowscanr)


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

### Include DEST sets
samps <- rbind(samps[set=="DrosRTEC"],
               samps[set=="DrosEU"],
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
  mutate(SNP_id = paste(chr, pos, sep = "_"))

snps.dt %<>%
  filter(SNP_id %in% colnames(dat_for_Analysis))

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R")]$variant.id)

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

###################
### Next part
### Calculate FST

#Generate outfile object
outfile = data.frame(
  samp1 = rep(NA, dim(comp_vector_samepops)[1]),
  samp2 = rep(NA, dim(comp_vector_samepops)[1]),
  FST = rep(NA, dim(comp_vector_samepops)[1])
)


for(i in 1:dim(comp_vector_samepops)[1]){
  
  print(i/dim(comp_vector_samepops)[1] * 100)
  
  samps_to_compare = c(comp_vector_samepops$samp1[i], comp_vector_samepops$samp2[i])
  
  pool_sizes = c(samps_YwY_EC$nFlies[which(samps_YwY_EC$sampleId == comp_vector_samepops$samp1[i])],
                 samps_YwY_EC$nFlies[which(samps_YwY_EC$sampleId == comp_vector_samepops$samp2[i])])
  
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare ),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
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
  
  
}## i   


left_join(comp_vector_samepops, outfile) -> Out_comp_vector_samepops


samps %<>%
  mutate(month_col = month(as.Date(collectionDate, 
                                   format = c("%m/%d/%Y"))))

samps %>%
  select(sampleId, month_col, year) -> samp_1_meta
names(samp_1_meta) = c("samp1", "month1", "year1")

samps %>%
  select(sampleId, month_col, year) -> samp_2_meta
names(samp_2_meta) =  c("samp2", "month2", "year2")

Out_comp_vector_samepops %<>%
  left_join(samp_1_meta) %>% 
  left_join(samp_2_meta) 

Out_comp_vector_samepops %<>%
  mutate(bin_date = ifelse(.$day_diff <= 200, "1.within", 
                           ifelse(.$day_diff >= 550, "3.Multi-Year", "2.Overwinter" ) ))

#Out_comp_vector_samepops %<>%
#	  .[which(.$pop1 %in%  cities_for_y_to_y),] 

head(Out_comp_vector_samepops)

### Save object
print("Save object")
save(Out_comp_vector_samepops, file = "Year_to_year_object.Rdata" )

### Save object
load("./Year_to_year_object.Rdata")

#--> object is Out_comp_vector_samepops
#---> also comp_vector

### Panel B
Out_comp_vector_samepops %>%
  .[which(.$pop1 == "Linvilla" & 
            .$pop2 == "Linvilla" |
            .$pop1 == "Charlottesville" & 
            .$pop2 == "Charlottesville"
  ),] %>%  
  ggplot(
    aes(
      x=day_diff,
      y=(FST/(1-FST)),
      fill=pop1,
      shape=bin_date
    ))  + 
  theme_bw() +
  geom_point(	alpha = 0.8,
              size = 2) +
  geom_smooth(method = "lm",
              color = "black") +
  xlim(0,1500) +
  xlab(expression(delta[Time])) +
  ylab(expression(F[ST]/1-F[ST])) +
  scale_shape_manual(values = c(21,22, 23)) +
  scale_fill_manual(values = c("thistle","grey73")) ->
  fst_Linv_Cville

ggsave(fst_Linv_Cville,
       file = "fst_Linv_Cville.pdf",
       width = 5,
       height = 3)

### Panel C

#selected_cities <-
#  c("Linvilla",
#    "Cross Plains",
#    "Charlottesville",
#    "Munich",
#    "Broggingen",
#    "Akaa")

Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  group_by(pop1, bin_date) %>% 
  summarise(FST_mean = mean(FST)) %>%
  dcast(pop1 ~ bin_date) ->
  mean_fst
  
Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%  
  left_join(mean_fst) %>% 
  ggplot(
    aes(
      x=fct_reorder(pop1, `1.within`),
      y=FST,
      color=bin_date
    )
  ) +
  geom_boxplot(width = 0.7) +
  xlab("Sampling Locale") +
  coord_flip() +
  theme_bw() +
  ylab(expression(italic(F)[ST])) -> 			
  fst_boxplot

ggsave(fst_boxplot,
       file = "fst_boxplot.pdf",
       width = 5,
       height = 4)

##########
Out_comp_vector_samepops %>%
  .[which(.$bin_date %in% 
            c("2.Overwinter", "1.within") ),] %>%
  dcast(pop1+samp1+samp2~bin_date, value.var = "FST") ->
  comp_vector_for_t

####

o_list = list()
  
comp_vector_for_t$pop1 %>% unique -> select_pop

for(i in 1:length(select_pop)){
  
tmp <- comp_vector_for_t %>%
        filter(pop1 == select_pop[i])

t.test(tmp$`1.within`[complete.cases(tmp$`1.within`)],
       tmp$`2.Overwinter`[complete.cases(tmp$`2.Overwinter`)]) ->
  tmp_a

o_list[[i]] = data.frame(pop=select_pop[i], P_val= tmp_a$p.value)

}

o_df = do.call(rbind, o_list)


