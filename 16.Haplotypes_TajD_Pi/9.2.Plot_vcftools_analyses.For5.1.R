### Explore pi and D in w5.1
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)
library(patchwork)
library(data.table)
library(SNPRelate)
library(SeqArray)

### GLM outliers

### Temp SNPs
models = c("temp.max;2;5.Cville" #,
           #"temp.ave;9;3.Europe_E",
           #"humidity.ave;4;2.North_America_W",
           #"humidity.ave;8;1.Europe_W"
)

#### base files
base <- "/project/berglandlab/alan/environmental_ombibus_global"
k=1

#all.wins.out = 
#foreach(k=1:length(models), .combine = "rbind")%do%{
#

file <- paste(base, models[k], paste(models[k],"glmRNP.Rdata", sep = ".") , sep = "/" )
print(file)

message(models[k])
out.glm <- get(load(file))

out.glm %<>%
  group_by(perm) %>%
  mutate(p_ltr_adj = p.adjust(p_lrt, method = "fdr"))

out.glm %>%
  filter(chr == "2L") %>%
  #filter(perm == 0) %>%
  filter(pos > 5.05e6-1e5, pos < 5.25e6+1e5) ->
  out.glm.w5.1

##
out.glm.w5.1 %>%
  mutate(perm.type = case_when(perm == 0 ~ "real",
                   perm != 0 ~ "perm" )) %>%
  group_by(perm.type, pos) %>%
  summarise(uci = quantile(rnp, 0.000001)) %>%
  dcast(pos ~ perm.type) %>% 
  mutate(sig = case_when(real < perm ~ "sig",
                         TRUE ~ "ns")) %>%
  filter(sig == "sig") -> glm_snps

### Mps 300 gff
gff.msp.300 = vroom("gff.msp.300.dat.txt", col_names = F, delim = "\t")
for(i in 1:dim(gff.msp.300)[1] ){
  tmp <- gff.msp.300[i,]
  
  tmp.str = str_split(tmp, pattern = "\\|")
  
  isoform = gsub('.+(isoform .).+','\\1',tmp.str[[9]][3])
  
  gff.msp.300$isoform[i] = isoform
  
}
gff.msp.300$isoform %>% unique()

gff.msp.300 %<>% 
  mutate(rank = case_when(isoform == "isoform B" ~ -0.001,
                          isoform == "isoform D" ~ -0.002,
                          isoform == "isoform E" ~ -0.003,
                          isoform == "isoform F" ~ -0.004,
                          isoform == "isoform G" ~ -0.005,
                          isoform == "isoform H" ~ -0.006,
                          isoform == "isoform I" ~ -0.007,
                          isoform == "isoform J" ~ -0.008,
                          isoform == "isoform K" ~ -0.009,
                          isoform == "isoform L" ~ -0.010,
                          isoform == "isoform M" ~ -0.011
  ))
#####


files <- system(" ls | grep 'win5.1' | grep 'pi' | grep -v 'pdf' ", intern = T)

pi.d.values <- foreach(i=1:length(files), .combine = "rbind")%do%{
  
  tmp = vroom(files[i])
  
  #if( grepl("Tajima", files[i])  ){
  #  tmp %>% melt(id = c("CHROM", "BIN_START", "N_SNPS") ) -> tmp.melt
  #  tmp.melt %<>%
  #    mutate(BIN_END = BIN_START+10000)
  #    names(tmp.melt)[3] = "N_VARIANTS"
  #} else
  if( grepl("pi", files[i])  ){
      tmp %>% reshape2::melt(id = c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS") ) -> tmp.melt
    }
  
  str_split(files[i], pattern = "\\.") -> metadat.tmp
  
  tmp.melt %<>%
    mutate(pop = metadat.tmp[[1]][4],
           kar = metadat.tmp[[1]][7])
  
  return(tmp.melt)
  
}

pi.d.values %<>% filter(BIN_START > 5.05e6-1e5, BIN_END < 5.25e6+1e5)
pi.d.values %>%
  filter(variable == "PI")  ->
  pi.2L.scales

##load FST
fst.cm <- vroom("inv.vs.std.cm.fst.windowed.weir.fst", na = c("", "NA", "-nan"))

fst.cm %>%
  filter(CHROM == "2L") ->
  fst.all.2L
fst.all.2L %<>% filter(BIN_START > 5.05e6-1e5 & BIN_END < 5.25e6+1e5)

fst.all.2L %>%
  arrange(-MEAN_FST)


#### Age


setwd("")

Dmel_name <- system("ls ../9.GEVA | grep 'Dmel' ", intern = T)

age.matrix <-
  foreach(i=Dmel_name 
          , .combine = "rbind"
  )%do%{
    
    message(i)
    
    tmp <- fread( 
      paste("../9.GEVA/",i, "/", paste(i, ".0.01.TMRCA.txt", sep =""), 
            sep = "" )
    )        
    tmp
  }

age.matrix %>%
  filter(V11 == "Dmel_ALL") %>%
  separate(V1, into = c("chr", "array_start", "array_end")) %>%
  mutate(SNP_id = paste(chr, V2, sep = "_") ) ->
  age.all


###### Plot
colors = c("purple", "#F8766D", "#00BA38", "#00BFC4", "#FF61CC")
### PT1
pi.2L.scales %>% 
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=value,
    color = pop,
    linetype = kar
  )) +
  scale_color_manual(values = colors) +
  ylab(expression(pi)) +
  geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  theme_classic() +
  xlab("Genomic Position (Mb)") +
  geom_line() ->
  pi.dat.5
ggsave(pi.dat.5, file = "pi.dat.5.pdf", w = 6, h =3)

####
fst.all.2L %>% 
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=MEAN_FST,
    #color = pop,
    #linetype = kar
  )) +
  geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  scale_color_manual(values = colors) +
  ylab(expression(paste("Wier and Cockerham ", F[ST] ) )) +
  theme_classic() +
  xlab("Genomic Position (Mb)") +
  geom_line() ->
  fst.dat.5
ggsave(fst.dat.5, file = "fst.dat.5.pdf", w = 6, h =3)
###

age.all %>%
  filter(V2 > 5.05e6-1e5 & V2 < 5.25e6+1e5) %>%
  ggplot(aes(
    x=V2,
    y=V8
  )) +
  geom_point() +
  geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) ->
  geva.5.1
ggsave(geva.5.1, file ="geva.5.1.pdf")  

##
ggplot() +
  #geom_point(data = glm_snps,
  #           aes(
  #             x=pos,
  #             y=0
  #           ), color = "red", shape = 5, size = 0.6) +
  geom_rect(data = gff.msp.300,
            aes(xmin=X4, xmax = X5,
                ymin = rank-0.0001, ymax = rank+0.0001), 
            alpha = 0.7, fill = "red") +
  geom_text(
    data = slice_head(group_by(gff.msp.300, isoform)),
    aes(x=5.07e6, y = rank,
    label = isoform
  ), size = 2.5 ) +
  theme_classic() +
  xlab("Genomic Position (Mb)") +
  geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  xlim(5.05e6-1e5,5.25e6+1e5) ->
  msp300.arc
  
ggsave(fst.dat.5/msp300.arc, file = "pi.fst.dat.5.pdf", w = 6, h =7)


###
### samps

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
glm_snps %>%
  filter(pos >= 5190001 & pos <= 5200000)

age.all %>%
  filter(V2 == 5192177)
#1 5192166 0.0033685604 3.347436e-03 sig -- sym
#2 5192177 0.0009666938 1.978267e-04 sig -- missense_variant**
#3 5192199 0.0001667910 4.685369e-05 sig -- sym
#4 5195810 0.0237046347 5.867123e-03 sig -- utr
#5 5196174 0.0062474470 1.601702e-03 sig -- sym
#6 5197284 0.0050540873 5.067140e-04 sig -- down
#7 5197312 0.0132465294 7.940832e-03 sig -- down
#8 5198138 0.0065250537 1.766558e-03 sig -- ups

getData(chr="2L", start=5192177, end=5192177) -> msp300.target

load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/nasa_power.weather.mod.Rdata")
names(weather.ave)[1] = "sampleId"

sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))

weather.ave %>%
  filter(mod == 2) %>%
  dplyr::select(sampleId, temp.max) -> tmax





msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  ggplot(aes(
    x=temp.max,
    y=af,
    color=as.factor(year)
  )) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black") ->
  msp300.exa

ggsave(msp300.exa, file = "msp300.exa.pdf")


 

