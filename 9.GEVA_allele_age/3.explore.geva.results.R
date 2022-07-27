## GEVA
## 

library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(foreach)
library(SeqArray)
#library(glmnet)
library(doMC)
registerDoMC(5)

Dmel_name <- system("ls | grep 'Dmel' ", intern = T)

age.matrix <-
foreach(i=Dmel_name 
        , .combine = "rbind"
        )%do%{
  
  message(i)
  
  tmp <- fread( 
    paste(i, "/", paste(i, ".0.01.TMRCA.txt", sep =""), sep = "" )
    )        
  
  tmp
  
}


### Load dat
### 
#geva.dat <- fread("./Dmel_GEVA.0.01.TMRCA.txt")

age.matrix %<>%
  separate(V1, into = c("chr", "array_start", "array_end")) %>%
  mutate(SNP_id = paste(chr, V2, sep = "_") )

###############
### windows ###
###############
# generate a master index for window analysis
### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(geva.dat, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L"
                        #,"2R", "3L", "3R"
                        ),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- age.matrix %>%
                    filter(chr == chr.i & V11 == "Dmel_ALL")
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$V2), to=max(tmp$V2)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$V2), to=max(tmp$V2)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
#####
#####
#####

wins %<>% filter(chr == "2L")
### start the summarization process
win.out.geva <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  win.tmp <- age.matrix %>%
              filter(chr == wins$chr[win.i]) %>%
              filter(V2 >= wins$start[win.i] & V2 <= wins$end[win.i] )
  
  win.tmp %>%
    group_by(V11) %>%
    summarise(mean.pos = mean(V2),
              med.age = median(V8),
              uci = quantile(V8, 0.975),
              lci = quantile(V8, 0.025)) ->
    data.tmp
  
  data.frame(wins$chr[win.i], data.tmp)
  
} ## close do




final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start, xmax = end,
                ymin = 0, ymax = 190000), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=filter(win.out.geva, V11 != "Dmel_ALL"),
    aes(
    x=mean.pos,
    y=(med.age*0.06666667),
    color = V11
  )) +
  #geom_ribbon() +
  geom_line() +
  ylab("TMRCA gens (x1e6)") +
  xlab("Genomic Position (Mb)") +
  ggtitle("TMRCA years") +
  theme_bw() +
  xlim(0,20.5*1e6) +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  geva.win

ggsave(geva.win, file ="geva.win.pdf", w = 7, h = 3)

win.out.geva %>% 
  filter(V11 == "Dmel_inv" & mean.pos > 6e6 & mean.pos <7e6) %>% 
  arrange(med.age) %>%
  mutate(age.est = med.age*0.06666667,
         uci.est = uci*0.06666667,
         lci.est = lci*0.06666667 )


####
age.matrix %>%
  filter(V11 == "Dmel_inv" & V2 > 6671911-1.0e5 & V2 < 6671911+1.0e5) %>% 
  ggplot(aes(
    x=V2/1e6,
    y=log10(V8*0.06666667),
  )) +
  geom_point(size = 0.9) +
  ggtitle("SNP-wise allele age") +
  xlab("Genomic Position (Mb)") +
  ylab(expression(Log[10](Years))) +
  geom_smooth(span = 1/10) ->
  sweep.age
ggsave(sweep.age, file ="sweep.age.pdf", w = 4, h = 3)

## find sweep target
## 
## 
## 
age.matrix %>%
  filter(V11 == "Dmel_inv" & V2 > 6597656-2e5 & V2 < 6597656+2e5) %>% 
  arrange(V8)

### find annotation
###
###
###
###
### samps
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")

### open GDS for common SNPs (PoolSNP)

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

annot.file <- seqGetData(genofile, "annotation/info/ANN")



##create get data function
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

### Implement function
### 
getData(chr="2L", start=6671911, end=6671911)

### get extra info
###load("/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata")
