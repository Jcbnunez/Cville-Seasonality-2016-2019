#### Calculate association of snps with inversion

library(tidyverse)
library(data.table)
library(vcfR)
library(foreach)
library(magrittr)

samps <- fread(
"/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Single_Individuals/TableS8.SVM.predictions.txt", 
header = T)

samps %<>%
mutate(INVGENO = case_when(
INV_STATUS == "STD" ~ 0,
INV_STATUS == "HET" ~ 1,
INV_STATUS == "INV" ~ 2))

###
###
###

obj.vcfR <- read.vcfR("CM.2L.maf05.mis85.thin100.inv.std.recode.vcf.gz")

geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information

SNP_names = data.frame(chromosome, position) %>%
mutate(SNP_id = paste(chromosome, position, sep = "_"))

obj.vcfR@gt %>% colnames %>% .[-1] -> sampleIds

###colnames(G)
### 9 here is the base code for NA
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

table(as.vector(G))
colnames(G)= sampleIds
rownames(G)= SNP_names$SNP_id

#####
#2L_2051609
#i=which(rownames(G) == "2L_2226674")
isEmpty <- function(x) {
    return(length(x)==0)
}


r2.invsnp = 
foreach(i=1:dim(G)[1],
.combine = "rbind",
.errorhandling = "remove")%do%{


G[i,] %>%
data.frame(geno.snp = .) %>%
mutate(sampleId = rownames(.)) %>%
left_join(samps) %>%
mutate(geno.joint = paste(geno.snp,INVGENO, sep = "|")) %>%
.[complete.cases(.),] -> tmp

### chrom number
Ngenos = dim(tmp)[1]*2

### allele frequencies
fA = sum(tmp$geno.snp)/Ngenos
fB = sum(tmp$INVGENO)/Ngenos

### Genotype freqeuncies

tmp$geno.joint %>% table %>% data.frame -> geno.tab
names(geno.tab)[1] = "genocall"


p00 = filter(geno.tab, genocall == "0|0")$Freq
p10 = filter(geno.tab, genocall == "1|0")$Freq
p22 = filter(geno.tab, genocall == "2|2")$Freq
p12 = filter(geno.tab, genocall == "1|2")$Freq
p02 = filter(geno.tab, genocall == "0|2")$Freq
p20 = filter(geno.tab, genocall == "2|0")$Freq

####
###
if(isEmpty(p00)){
p00 = 0
}
if(isEmpty(p10)){
p10 = 0
}
if(isEmpty(p22)){
p22 = 0
}
if(isEmpty(p12)){
p12 = 0
}
if(isEmpty(p02)){
p02 = 0
}
if(isEmpty(p20)){
p20 = 0
}


####
pAB = p00*2 + p10
pab = p22*2 + p12
pAb = p02*2 + p12
paB = p20*2 + p10

###
data.frame(
SNPid = rownames(G)[i],
p00,p10,p22,p12,p02,p20,
fAB=pAB/Ngenos,
fab=pab/Ngenos,
fAb=pAb/Ngenos,
faB=paB/Ngenos,
fA = fA,
fB = fB
) %>%
mutate(D = fAB*fab-fAb*faB) %>%
mutate(Dmax = fA*(1-fB)) %>%
mutate(Dprime = D/Dmax) %>%
mutate(r2 = (D^2)/(fA*(1-fA)*fB*(1-fB))) -> o

#see https://doi.org/10.1159/000504171

message(paste(i,dim(G)[1],o$r2, sep = "|"))

return(o)
}

save(r2.invsnp, file = "r2.invsnp.CM.Rdata")

r2.invsnp %>%
as.data.frame() %>%
separate(SNPid, into = c("chr","pos"), sep = "_") %>%
ggplot(aes(
x=as.numeric(pos),
y=r2
)) +
geom_point() +
geom_smooth() ->
r2.plot
ggsave(r2.plot, file = "r2.plot.pdf")

