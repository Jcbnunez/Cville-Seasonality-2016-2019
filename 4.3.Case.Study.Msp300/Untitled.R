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

####

seqVCF2GDS("./Simcline.Dsim_Scf_2L.4996892.recode.vcf", 
           "./Simcline.Dsim_Scf_2L.4996892.recode.gds", 
           storage.option="ZIP_RA")

####
### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("./Simcline.Dsim_Scf_2L.4996892.recode.gds")

snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"),
                     sra = seqGetData(genofile, "sample.id"),
                     Dosage = seqGetData(genofile, "$dosage")
                     )

### add metadata
metadat <- vroom("concatenated.csv")

### Join

left_join(snp.dt, metadat) %>%
  .[complete.cases(.$Dosage.V1),] %>%
  group_by(city, country) %>%
  summarise(lat = mean(lat), long = mean(long), Allele.count= sum(Dosage.V1), n.ind = n() ) %>%
  mutate(A = Allele.count, B = (n.ind*2)-(Allele.count) ) %>%
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


