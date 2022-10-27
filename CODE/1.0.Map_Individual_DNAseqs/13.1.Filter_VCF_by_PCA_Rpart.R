###########################################################################
###########################################################################
# Load R libraries
###########################################################################
###########################################################################

library(vcfR)
library(adegenet)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(magrittr)
library(reshape2)
library(patchwork)

###########################################################################
###########################################################################
# Load data
###########################################################################
###########################################################################

filter_vcf = "/scratch/yey2sn/phasing_droso_data/Whatshap/PCA_filter/PrePCA.mis80_minQ100_minDP10_maf5_thin1k.3L.recode.vcf.gz"

#Load VCF
working_VCF = read.vcfR(filter_vcf)

#Transform VCF to genlight
Genlight <- vcfR2genlight(working_VCF)

#load Sample names
sample_names <- read.table("/scratch/yey2sn/phasing_droso_data/Whatshap/PCA_filter/PrePCA.mis80_minQ100_minDP10_maf5_thin1k.3L.sample_names.txt")

#Generate some metadata
sample_names %<>% separate(V1, into = c("Plate","Well","Pop"), remove = F)
sample_names$Plate = gsub("[0-9]","", sample_names$Plate)
sample_names$Pop[which(sample_names$Plate == "CM")] = "TYS"

#Include sample/pop names into the genlight
Genlight@pop=as.factor(sample_names$Pop)

#Build PCA

natural_samps = sample_names %>%
					.[which(.$Pop %in% c("TYS", "CMspring","CMfall")),]


write.table(natural_samps$V1,
			file = "natural_pops_sampleid.txt",
			sep = "\t",
			quote = F,
			row.names = F,
            col.names = F
            )


tab(Genlight) %>% 
  .[which(rownames(.) %in% natural_samps$V1), ] %>%
	PCA(graph = F,scale.unit = F) -> test_PCA

test_PCA %>% fviz_pca_ind(habillage = as.factor(natural_samps$Pop),
geom = c("point"), 
addEllipses = T, 
repel = T  ) -> test_PCA_g
ggsave("test_PCA_g.pdf", test_PCA_g)


test_PCA %>% fviz_pca_ind(habillage = as.factor(natural_samps$Pop),
geom = c("point"), 
addEllipses = T,
axes = c(2,3), 
repel = T  ) -> test_PCA_g23
ggsave("test_PCA_g23.pdf", test_PCA_g23)

