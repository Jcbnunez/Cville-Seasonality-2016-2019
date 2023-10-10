### Create final object
library(tidyverse)

pheno_dat <- get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.Nunez_et_al_Supergene_paper/TTsnps.cleaned"))

geno_dat <- get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Revision_Best_Models/temp.max;2;5.Cville.v2.glmRNP.Rdata"))

### ==> Get columns of interest
#chr	chromosome
#pos	position
#AIC	AIC of the model at given SNP
#variable	environmental variable from NASA-POWER
#mod_id	Model Id
#p_lrt	LRT p value
#b_temp	beta of environemtnal variable
#se_temp	SE of environemtnal variable
#Cluster	Population, i.e., VA, EUE, EUW
#cm_mb	Recombiantion
#invName	Inversion Name
#rnp	Ranked Normalized P-value
#Perm.rnp.0.01.quant	0.01 quantile of RNP in 100 permutations
#time_window	window of time of the model
#SNP id	SNP id
N.phenos	Phenotype number
Description	Pheno description

####
geno_dat %>%
group_by(perm==0,chr, pos) %>%
summarize(Perm.rnp.0.01.quant=quantile(rnp, 0.01)) %>%
filter(`perm == 0` == TRUE) %>%
mutate(SNP_id = paste(chr, pos, sep = "_")) -> rnp.summ
####
rnp.summ %>%
dplyr::select(SNP_id, Perm.rnp.0.01.quant) -> part2.obj
part2.obj[,-c(1:2)] -> part2.obj

####
geno_dat %>%
filter(perm == 0) %>%
dplyr::select(chr, pos, AIC, variable, mod_id=mod, p_lrt=p_lrt.x,
b_temp, se_temp, cluster, cm_mb, invName=inv, rnp ) %>%
mutate(time_window = "0-15d", 
SNP_id = paste(chr, pos, sep = "_") ) -> part1.obj
####
pheno_dat %>% 
mutate(SNP_id = paste(chr,pos, sep = "_"))%>%
dplyr::select(SNP_id, N.phenos=pheno.n, 
Description= pheno.list) -> part3.obj

####
part1.obj %>%
left_join(part2.obj) %>% 
left_join(part3.obj)-> DataS1



