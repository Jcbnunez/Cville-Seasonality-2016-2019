### Panels 3B and 3C
### 
rm(list = ls())


library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)
library(SNPRelate)
library(SeqArray)
library(vroom)
library(fastglm)

### Panel 3B
### 

#save(sub_pi_d_parsed.plot, inv.dt, final.windows.pos, file = "dat.for.3b.Rdata")
load("./dat.for.3b.Rdata")

ggplot() + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="solid") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="solid") +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = -0, ymax = 0.01), 
            alpha = 0.7, fill = "gold") +
  geom_line(
    data=sub_pi_d_parsed.plot,
    aes(
      x=BIN_START/1e6,
      y=value,
      color = type),
    alpha = 0.9) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,20.5) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") ->
  pi_d_plot_all

ggsave(pi_d_plot_all, file = "pi_d_plot_all.pdf", w = 7, h = 3.5)

### Panel 3C
### 

fst <- vroom("inv.vs.std.cm.fst.windowed.weir.fst")

fst %>% 
  filter(CHROM == "2L") -> fst.2l.dat
  
fst.2l.dat %>%
  arrange(-MEAN_FST)
  
  final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start, xmax = end,
                ymin = 0, ymax = 0.7), 
            alpha = 0.7, fill = "gold") +
  geom_line(data= fst.2l.dat, aes(
    x=((BIN_START+BIN_END)/2),
    y=MEAN_FST,
    #color = pop,
    #linetype = kar
  ),  color = "purple", alpha = 0.9) +
  #geom_vline(xintercept = (5170001+5180000 )/2, color = "blue" ) +
  #geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  scale_color_manual(values = colors) +
  ylab(expression(paste(F[ST]))) +
  theme_bw() +
  xlab("Genomic Position (Mb)") +
  xlim(0,20.5*1e6) +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  fst.dat.2l
ggsave(fst.dat.2l, file = "fst.dat.2l.pdf", w = 7, h = 3.5)

### Panel zoom on msp300
### 


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
fst.cm %>%
  filter(CHROM != "X") %>%
  ggplot(aes(
    x=((BIN_START+BIN_END)/2)/1e6,
    y=MEAN_FST,
  )) +
  geom_line() +
  ylab(expression(F[ST])) +
  xlab("Genomic Position") +
  facet_wrap(~CHROM, scales = "free_x") ->
  all.fasts
ggsave(all.fasts, file = "all.fasts.pdf", w= 7, h = 5)
ggsave(all.fasts, file = "all.fasts.png", w= 7, h = 5)

fst.cm %>%
  filter(CHROM == "2L") ->
  fst.all.2L

fst.all.2L %>%
  filter(BIN_START > 2225744 & BIN_END < 13154180) %>%
  summarise(Mean = mean(.$MEAN_FST, na.rm = T),
            SD = sd(.$MEAN_FST, na.rm = T))

fst.all.2L %>%
  arrange(-MEAN_FST) %>% head


fst.all.2L.w5.1 <- fst.all.2L %>% filter(BIN_START > 5.05e6-1e5 & BIN_END < 5.25e6+1e5)

fst.all.2L %>% 
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=MEAN_FST,
    #color = pop,
    #linetype = kar
  )) +
  geom_vline(xintercept = (5190001+5200000 )/2, color = "red" ) +
  scale_color_manual(values = colors) +
  ylab(expression(paste("Wier and Cockerham ", F[ST] ) )) +
  theme_classic() +
  xlab("Genomic Position (Mb)") +
  geom_line() ->
  fst.dat.5
ggsave(fst.dat.5, file = "fst.dat.5.pdf", w = 6, h =3)

ggplot() +
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


####
####

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
  # chr="2L"; start=5192177; end=5192177
  #chr="2L"; start=5192177; end=5192177
  
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


### get data
getData(chr="2L", start=5192177, end=5192177) -> msp300.target

msp300.target %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC")) %>%
  summarise(mean.af = mean(af, na.rm = T))

msp300.target %>%
  filter(set %in% c("CvilleSet", "DrosEU", "DrosRTEC")) %>%
  group_by(continent) %>%
  summarise(mean.af = mean(af, na.rm = T))

msp300.target %>%
  filter(sampleId == "SIM") 

#p.Gly10912Val

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
    y=af)
  ) + 
  theme_bw() +
  geom_smooth(method = "lm", fill = "grey", alpha = 0.3, color = "black") +
  geom_point(aes(shape=as.factor(year)), size = 3) ->
  msp300.exa

ggsave(msp300.exa, file = "msp300.exa.pdf", w = 4, h =3)


msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  ggplot(aes(
    x=as.Date(collectionDate, format = "%m/%d/%Y"),
    y=af,
    color = as.factor(year))
  ) + 
  theme_bw() +
  geom_point() +
  geom_smooth(method = "lm") ->
  #geom_smooth(method = "lm", fill = "grey", alpha = 0.3, color = "black") +
  #geom_point(aes(shape=as.factor(year)), size = 3) ->
  msp300.exa.time

ggsave(msp300.exa.time, file = "msp300.exa.time.pdf", w = 4, h =3)



#1 5192166 0.0033685604 3.347436e-03 sig -- sym
#2 5192177 0.0009666938 1.978267e-04 sig -- missense_variant**
#3 5192199 0.0001667910 4.685369e-05 sig -- sym
#4 5195810 0.0237046347 5.867123e-03 sig -- utr
#5 5196174 0.0062474470 1.601702e-03 sig -- sym
#6 5197284 0.0050540873 5.067140e-04 sig -- down
#7 5197312 0.0132465294 7.940832e-03 sig -- down
#8 5198138 0.0065250537 1.766558e-03 sig -- ups
### model

msp300.target %>%
  left_join(tmax) %>%
  filter(locality == "VA_ch") %>%
  filter(year %in% 2016:2018) %>%
  .[complete.cases(af),] -> dat.msp300.5192177

dat <- dat.msp300.5192177
glm.method <- 0
X.year.var <- model.matrix(~as.factor(year)+temp.max, dat)

t1.year.var <- fastglm(x=X.year.var, y=dat$af_nEff,
                       family=binomial(), weights=dat$nEff, 
                       method=glm.method)
summary(t1.year.var)
####
X.each.year.var <- model.matrix(~temp.max, dat)

dat <- filter(dat.msp300.5192177, year == 2016)
t2016.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                       family=binomial(), weights=dat$nEff, 
                       method=glm.method)
summary(t2016.year.var)
###
dat <- filter(dat.msp300.5192177, year == 2017)
X.each.year.var <- model.matrix(~temp.max, dat)
t2017.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                          family=binomial(), weights=dat$nEff, 
                          method=glm.method)
summary(t2017.year.var)
###
dat <- filter(dat.msp300.5192177, year == 2018)
X.each.year.var <- model.matrix(~temp.max, dat)
t2018.year.var <- fastglm(x=X.each.year.var, y=dat$af_nEff,
                          family=binomial(), weights=dat$nEff, 
                          method=glm.method)
summary(t2018.year.var)





