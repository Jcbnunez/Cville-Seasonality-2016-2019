####### Plot haplotypes
####### 
####### 
rm(list = ls())
### libraries
library(data.table)
library(SeqArray)
library(tidyverse)
library(car)
library(foreach)

### load meta-data file
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/project/berglandlab/DEST_Charlottesville_TYS/DEST_pipeline_output/dest.all.PoolSNP.001.50.10Mar2021.ann.gds")

### common SNP.dt
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))
snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, snp.dt$id)

### function
getData <- function(chr="2L", start=14617051, end=14617051) {
  # chr="2L"; start=14617051; end=14617051
  
  ### filter to target
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)
  
  ### get annotations
  #message("Annotations")
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
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  afi <- merge(af, snp.dt1.an, by="variant.id")
  afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis[,-c("n"), with=F]
}

### test
data <- getData()


### Regions
load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")
head(final_in2Lt_markers)
left_break = getData(chr="2L", start=2051609 , end=3096574)
right_break = getData(chr="2L", start=11584064 , end=13204668)


###Generate windows
#load("/scratch/yey2sn/Overwintering_ms/4.GML_plots/PEAKS_for_ANALYSIS.Rdata")
#PEAKS_for_ANALYSIS

#starts= c(5155762, 6255762, 9505762)
#ends= c(5255762, 6355762, 9605762)


#### load windows
#haplo_windows <- "/scratch/yey2sn/Overwintering_ms/4.2.env.omibus.mods/final.windows.pos.Rdata"
#load(haplo_windows)
#final.windows.pos

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_9.6" ),
             start = c(3000026, 4650065, 5050026, 6100321, 9500026),
             end = c(3150026, 4799922,  5250026, 6224905, 9650026)) 


rbind(
mutate(left_break, win = "left"),
mutate(right_break, win = "right"),
getData(chr="2L", start=3000026, end=3150026) %>% mutate(win="win_3.1"),
getData(chr="2L", start=4650065, end=4799922) %>% mutate(win="win_4.7"),
getData(chr="2L", start=5050026, end=5250026) %>% mutate(win="win_5.2"),
getData(chr="2L", start=6100321, end=6224905) %>% mutate(win="win_6.1"),
getData(chr="2L", start=9500026, end=9650026) %>% mutate(win="win_9.6")) ->
AF_wins_alldat


#######

### MAKE VA DATAFRAME
AF_wins_alldat %>%
  as.data.frame() %>%
  filter(locality == "VA_ch") %>% 
  group_by(pos, win) %>%
  summarize(AF_mean = mean(af_nEff, na.rm = T)) ->
  AF_dat_summ

### MAKE SIM DATAFRAME
AF_wins_alldat %>%
  as.data.frame() %>%
  filter(sampleId == "SIM") %>% 
  dplyr::select(sampleId, variant.id, chr, pos, af, col, gene) ->
  SIM_AF

AF_dat_summ %>%
  #.[complete.cases(.$AF_mean),] %>% 
  #filter(AF_mean > 0.01) %>% 
  mutate(SNP_id = paste("2L", pos, "SNP", sep = "_")) ->
  AF_dat_summ_id

save(AF_dat_summ_id, SIM_AF, file = "AF_dat_summ_id.dat.Rdata")

write.table(AF_dat_summ_id$SNP_id, file = "Loci_for_haplotype_plotting_windows.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


