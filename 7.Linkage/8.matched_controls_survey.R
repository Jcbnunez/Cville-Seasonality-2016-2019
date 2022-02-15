### Script 8 of the LD analysis.
### Matched control analyses
### The point of this script is to obtain a number of "matched controls for the LD+GLM outliers"
### and investigate their patterns of LD.
### Do the LD patterns that we obsrse.. are they expected by chance?
### 

### libraries
library(data.table)
library(SeqArray)
library(tidyverse)
library(foreach)
library(car)
library(DescTools)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

###
### Part 1. Obtain SNPs of interest
### load meta-data file
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

### bring in temperature data
load("/project/berglandlab/DEST_Charlottesville_TYS/weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"

### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)

### common SNP.dt
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))
snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, snp.dt$id)


### This is a core function which extract SNPs of interest
### function
getData <- function(chr="2L", start=14617051, end=14617051, annotate_priority=FALSE) {
  # chr="2L"; start=14617051; end=14617051
  
  if(abs(start-end) != 0 &  annotate_priority==TRUE ){
    stop(' NOTICE: annotate_priority=TRUE only works for single SNPs! ')
    }
  
  ### filter to target
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)
  
  ### get annotations --- UPDATED by JCBN on Feb 9, 2022
  ### Now the script pick the "unique" SNP annotation based on priority scoring
  if(annotate_priority==TRUE){
  message("Annotations")

  tmp_annot <-  data.frame(annot = seqGetData(genofile, "annotation/info/ANN")[[2]])
  
  #extract annotation and parition into categories
  tmp_annot %>%
    separate(annot,
             into = c(
               "Allele",
               "Annotation",
               "Putative_impact",
               "Gene_Name",
               "Gene_ID",
               "Feature_type",
               "Feature_ID",
               "Transcript_biotype",
               "Rank",
               "HGVS.c",
               "HGVS.p",
               "cDNA_position",
               "CDS_position",
               "Protein_position",
               "Distance_to_feature",
               "Notes"),
             sep = "\\|",
             convert = F
    ) %>%
    dplyr::select(Allele, Annotation, Putative_impact, Gene_Name, Gene_ID) %>%
    mutate(priority = case_when(Putative_impact == "MODIFIER" ~ 0,
                                Putative_impact == "LOW" ~ 1,
                                Putative_impact == "MODERATE" ~ 2,
                                Putative_impact == "HIGH" ~ 3
                                )) %>%
    slice_max(priority, with_ties = F)  -> tmp_annot_sep


  ## include annotaions as well as position and chromosome
  snp.dt1.an <- data.frame(chr=chr, 
                             pos=start:end,
                             variant.id = snp.dt[J(snp.tmp), nomatch=0]$id,
                             tmp_annot_sep
  )
  
  } #close annotations
  
  ##### Classic annotations
  if(annotate_priority==FALSE){
    ### get annotations
    message("Annotations")
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data
    
    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))
    
    # Extract data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
    
    # Collapse additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                          list(variant.id=id)]
    
    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]
  } ## close annotate classic
  
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp)[,1],   #$data
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  
  if(annotate_priority==TRUE){
  afi <- merge(af, snp.dt1.an, by="variant.id")
  afi <- merge(afi, snp.dt[,-c(1:2)], by.x="variant.id", by.y=c("id"))
  } 
  
  if(annotate_priority==FALSE){
    afi <- merge(af, snp.dt1.an, by="variant.id")
    afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")
    } 

  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by=c("sampleId") )
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### Add temperature
  message("add temperature")
  afis = right_join(afis, weather.ave)
  ### return
  return(afis)
}

### test
### data <- getData()
### 

### Part 2. Bring in the SNP targets
GLM_LD_outliers <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers.txt")
## this is generated in step 5 of the LD scripts ....

### Part 3. extract position
AF_GLM_annots <- foreach(i=1:dim(GLM_LD_outliers)[1],
                         .combine = "rbind",
                         .errorhandling = "remove")%dopar%{
  getData(chr=GLM_LD_outliers$chr[i], 
          start=GLM_LD_outliers$pos[i], 
          end=GLM_LD_outliers$pos[i],
          annotate=TRUE)
}

### Save Object
save(AF_GLM_annots, file = "GLM_LD_outlier_AFs_annot.Rdata")
load("./GLM_LD_outlier_AFs_annot.Rdata")

AF_GLM_annots %>%
  filter(!is.na(af)) %>%
  mutate(het = 2*af_nEff*(1-af_nEff)) %>%
  group_by(variant.id, chr, pos, Annotation) %>%
  summarise(mean_het = mean(het)) %>%
  mutate(het_bin = RoundTo(mean_het, 0.1, "floor")) ->
  GLM_LD_metadata_object

### Generate a controls
### GET inversion data

inversion_SNPs = getData(chr="2L", 
        start=2225744, 
        end=13154180,
        annotate_priority=FALSE)

inversion_SNPs %>%
  filter(!is.na(af)) %>%
  mutate(het = 2*af_nEff*(1-af_nEff)) %>%
  group_by(variant.id, chr, pos, col) %>%
  summarise(mean_het = mean(het)) %>%
  mutate(het_bin = RoundTo(mean_het, 0.1, "floor")) ->
  inversion_SNPs_het


## this makes a general object with heterozygocity and annotation 
save(inversion_SNPs_het, file = "./inversion_SNPS_hetbins_annot.Rdata")
load("./inversion_SNPS_hetbins_annot.Rdata")


## extract matched controls
inversion_SNPs_het %>%
  filter(!pos %in% GLM_LD_metadata_object$pos) ->
  inversion_SNPs_het_sans_GLMsLDs


######
matched_controls <- foreach(i=1:dim(GLM_LD_metadata_object)[1],
                            .errorhandling = "remove",
                            .combine = "rbind" )%dopar%{
  
  message(i)
                              
  pos_m =   GLM_LD_metadata_object$chr[i]
  chr_m =   GLM_LD_metadata_object$pos[i]
  Annot_m = GLM_LD_metadata_object$Annotation[i]
  Het_bin_m = GLM_LD_metadata_object$het_bin[i]
    
  inversion_SNPs_het_sans_GLMsLDs %>%
    filter(het_bin ==  Het_bin_m,
           col == Annot_m) %>%
    as.data.frame() %>%
    slice_sample(n = 200, replace = FALSE) %>%
    mutate(matched_to_chr = chr_m,
           matched_to_pos = pos_m,
           delta_bp = abs(matched_to_chr-pos)) %>%
    filter(delta_bp > 1e6 & delta_bp < 10e6) %>%
    slice_sample(n = 100, replace = FALSE) 
  
}


save(matched_controls, file = "macthed_controls_interm.file.Rdata")





inversion_SNPs_het_sans_GLMsLDs$col %>% table






############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
### Part 3. Some explaration
### 

##bring in high didtance pairs
vars_of_interest = fread("/scratch/yey2sn/Overwintering_ms/7.LD/0.7_GLM_LD_outliers.txt")
vars_of_interest %>% filter(is.na(inversion_marker)) %>% .$pos -> pos_of_int


AF_GLM_annots %>%
  filter(locality %in% c("VA_ch"),
         pos %in% pos_of_int,
         #Gene_Name %in% genes_of_int,
         !is.na(af),
         Putative_impact %in% c(#"LOW",
                                "MODERATE", 
                                "HIGH")
                             ) %>%
  mutate(het = 2*af_nEff*(1-af_nEff)) %>%
  mutate(temp_round = RoundTo(aveTemp/10, 2, "floor") ) %>%
  mutate(month_summary = month(as.Date(collectionDate, format = "%m/%d/%Y")) ) -> in_dat

in_dat %>% 
  group_by(year, temp_round, variant.id) %>%
  summarise(mean_het = mean(het),
            mean_af = mean(af_nEff)) %>%
  mutate(type = "temp") ->  temp_summaries
names(temp_summaries)[2] = "summary"


in_dat %>% 
  group_by(year, month_summary, variant.id) %>%
  summarise(mean_het = mean(het),
            mean_af = mean(af_nEff)) %>%
  mutate(type = "month") -> month_summaries
names(month_summaries)[2] = "summary"

#rbind(temp_summaries, month_summaries) -> in_dat_r

month_summaries %>%
  ggplot(aes(
    x=summary,
    y=logit(mean_af),
    color = as.factor(variant.id),
  )) +
  geom_line(alpha = 0.9) +
  scale_color_manual(values = rep("grey", dim(month_summaries)[1])) +
  theme(legend.position = "none")  +
  facet_grid(year~type, space = "free") +
  ggtitle("Heterozygocity across GLM+LD outliers in 2lt [only missense variants consider]") ->
  test111

ggsave(test111, file = "test111.png",
       width = 6,
       height = 7)

####
temp_summaries %>%
  ggplot(aes(
    x=summary,
    y=logit(mean_af),
    color = as.factor(variant.id),
  )) +
  geom_line(alpha = 0.9) +
  scale_color_manual(values = rep("grey", dim(month_summaries)[1])) +
  theme(legend.position = "none")  +
  facet_grid(year~type, space = "free") +
  ggtitle("Heterozygocity across GLM+LD outliers in 2lt [only missense variants consider]") ->
  test111

ggsave(test111, file = "test111.png",
       width = 6,
       height = 7)

