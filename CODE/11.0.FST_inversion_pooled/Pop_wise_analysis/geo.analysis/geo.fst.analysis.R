##### Comparison of spatial FST
##### 

setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/11.0.FST_inversion_pooled/Pop_wise_analysis/geo.analysis")

load("FST.geographical.Rdata")


####
fst.dat.geo.EC %>%
  filter(Continental_clusters == "1.Europe_W",
         SNP.set == "glm.snps") -> EUW.glm
fst.dat.geo.EC %>%
  filter(Continental_clusters == "1.Europe_W",
         SNP.set == "macthed.controls.Inv") -> EUW.cI
fst.dat.geo.EC %>%
  filter(Continental_clusters == "1.Europe_W",
         SNP.set == "macthed.controls.noInv") -> EUW.cnoI


kruskal.test(EUW.glm$FST, EUW.cI$FST )
kruskal.test(EUW.glm$FST, EUW.cnoI$FST )

####
fst.dat.geo.EC %>%
  filter(Continental_clusters == "3.Europe_E",
         SNP.set == "glm.snps") -> EUE.glm
fst.dat.geo.EC %>%
  filter(Continental_clusters == "3.Europe_E",
         SNP.set == "macthed.controls.Inv") -> EUE.cI
fst.dat.geo.EC %>%
  filter(Continental_clusters == "3.Europe_E",
         SNP.set == "macthed.controls.noInv") -> EUE.cnoI


kruskal.test(EUE.glm$FST, EUE.cI$FST )
kruskal.test(EUE.glm$FST, EUE.cnoI$FST )
