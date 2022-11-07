#### Explore Courntey;s VCF
rm(list = ls())

library(foreach)
library(SeqArray)
library(gdsfmt)
library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(vroom)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)
library(vcfR)
library(scatterpie)

####
####
#module load bcftools
#module load vcftools
#
#bcftools query -l Simcline_final_2021.vcf.gz
#
#module load 
#vcftools \
#--gzvcf Simcline_final_2021.vcf.gz \
#--recode \
#--recode-INFO-all \
#--out Simcline.Dsim_Scf_2L.4996892  \
#--chr Dsim_Scf_2L \
#--from-bp 4996892 \
#--to-bp 4996892

#seqVCF2GDS("./Simcline.Dsim_Scf_2L.4996892.recode.vcf", 
#           "./Simcline.Dsim_Scf_2L.4996892.recode.gds", 
#           storage.option="ZIP_RA")

####
### open GDS for common SNPs (PoolSNP)
genofile.sim <- seqOpen("./Simcline.Dsim_Scf_2L.4996892.recode.gds")

snp.dt.sim <- data.table(chr=seqGetData(genofile.sim, "chromosome"),
                     pos=seqGetData(genofile.sim, "position"),
                     nAlleles=seqGetData(genofile.sim, "$num_allele"),
                     id=seqGetData(genofile.sim, "variant.id"),
                     sra = seqGetData(genofile.sim, "sample.id"),
                     Dosage = seqGetData(genofile.sim, "$dosage")
                     )

### add metadata
metadat <- vroom("concatenated.csv")

### Join

left_join(snp.dt.sim, metadat) %>%
  .[complete.cases(.$Dosage.V1),] %>%
  group_by(city, country) %>%
  summarise(lat = mean(lat), long = mean(long), Allele.count= sum(Dosage.V1), n.ind = n() ) %>%
  mutate(B = Allele.count, A = (n.ind*2)-(Allele.count) ) %>%
  mutate(AF = Allele.count/(2*n.ind) ) -> dat.for.plot.sim
  
####

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "antiquewhite", alpha = 0.8) +
  coord_sf(xlim = c(-125.15, 60.00), ylim = c(-38.00, 65.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_scatterpie(aes(x=long, y=lat,  r=2.5, group =country), data=dat.for.plot.sim,
                  cols=LETTERS[1:2], 
                  color="black") +
  xlab("Lon") + 
  ylab("Lat") + 
  ggtitle("D. Simulans Msp300 mutation") + 
  theme(legend.position = "none") -> sim.world
  
  ggsave(sim.world, file ="sim.world.pdf")

  ####   ####   ####   ####   ####   #### 
  ####     ####   ####   ####   ####   #### 
  ####       ####   ####   ####   #### 
  ####         ####   ####   ####   ####   #### 
 
  
   #### Melanogaster data
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
  
msp.300.dat = getData(chr="2L", start=5192177, end=5192177)

save(msp.300.dat, world, file= "msp.300.wolrd.plot.dat.Rdata")

msp.300.dat %>%
  filter(set %in% c("CvilleSet", "DrosEu", "DrosRTEC", "dgn" )) %>%
  group_by(continent) %>%
  summarise(af = mean(af))

msp.300.dat %>%
  filter(set %in% c("CvilleSet", "DrosEu", "DrosRTEC", "dgn" )) %>%
  separate(locality,
           into= c("region", "cit")) %>%
  .[complete.cases(af),] %>%
  group_by(region) %>%
  summarise(lat = mean(lat), long = mean(long), af = mean(af)) %>%
  mutate(B= round(af*100), A=round((1-af)*100)) ->
  msp.300.dat.plot

ggplot(data = world) +
  geom_sf(fill= "antiquewhite", alpha = 0.8) +
  coord_sf(xlim = c(-125.15, 60.00), ylim = c(-38.00, 65.00), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + 
  geom_scatterpie(aes(x=long, y=lat,  r=2.5), data=msp.300.dat.plot,
                  cols=LETTERS[1:2], 
                  color="black") +
  xlab("Lon") + 
  ylab("Lat") + 
  ggtitle("D. melanogaster Msp300 mutation") + 
  theme(legend.position = "none") -> mel.world

ggsave(mel.world, file ="mel.world.pdf")

