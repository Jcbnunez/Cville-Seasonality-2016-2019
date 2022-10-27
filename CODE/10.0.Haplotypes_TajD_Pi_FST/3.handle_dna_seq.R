### Analyse DNA sequencies
### 


library(Biostrings)
library(tidyverse)
library(magrittr)
library(data.table)

#record the fasta file address
#
type="Homozyg_only"
file_loc <- "/scratch/yey2sn/Overwintering_ms/11.Haplotypes/STD_INV.CM.joint.fasta"

#type="Het_only"
#file_loc <- "/scratch/yey2sn/Overwintering_ms/11.Haplotypes/HET.CM.joint.fasta"


#read in the fasta file
inv2lt_seq <-
  readDNAMultipleAlignment(filepath = file_loc,
                           format="fasta")

#load mask

#the begining of inv2Lt
pos_adj_scalar = 2051600

mask_info.all <- fread("./Retain_loci_metadat.txt")
mask_info.all %<>%
  mutate(adj_pos = V2-pos_adj_scalar+1)

mask_info.all %>%
  .$V3 %>% table


marker="GlmInv"
mask_info_glm = mask_info.all %>%
  filter(V3 %in% c("glm", "inv"))

mask_info_glm$adj_pos %>% table %>% table

colmask(inv2lt_seq, invert=TRUE) <- sort(mask_info_glm$adj_pos)

DNAStr_glm = as(inv2lt_seq, "DNAStringSet")
DNAStr_glm
writeXStringSet(DNAStr_glm,file=paste(type,marker,"fa", sep = "."))


### only get inverted sites
homozyg_samps <- fread("./Homozyg_samps_metadat.txt", header=F)
names(homozyg_samps) = c("samp","inv_st")
het_samps <- fread("./Het_samps_metadat.txt", header=F)
names(het_samps) = c("samp","inv_st")
rbind(homozyg_samps, het_samps) -> metadata_tab

metadata_tab %>%
  filter(inv_st == "Inv") -> samps_inv

karyo = "only_Invs"
only_Inv_DNA = inv2lt_seq
rowmask(only_Inv_DNA, invert=TRUE) <- grep(paste(samps_inv$samp, collapse="|"), rownames(inv2lt_seq))
DNAStr_only_Inv_DNA = as(only_Inv_DNA, "DNAStringSet")
DNAStr_only_Inv_DNA 
writeXStringSet(DNAStr_only_Inv_DNA,file=paste(type,marker,karyo,"fa", sep = "."))





#####individual groups
##marker="Glm"
##mask_info_glm = mask_info.all %>%
##  filter(V3 %in% c("glm"))
##glm_seq = inv2lt_seq
##colmask(glm_seq, invert=TRUE) <- mask_info_glm$adj_pos
##DNAStr_glm = as(glm_seq, "DNAStringSet")
##writeXStringSet(DNAStr_glm,file=paste(type,marker,"fa", sep = ".") )
##
##marker="Inv_left"
##mask_info_inv_l = mask_info.all %>%
##  filter(V3 %in% c("inv"),
##         V2 < 5e6)
##inv_l_seq = inv2lt_seq
##colmask(inv_l_seq, invert=TRUE) <- mask_info_inv_l$adj_pos
##DNAStr_inv_l = as(inv_l_seq, "DNAStringSet")
##writeXStringSet(DNAStr_inv_l,file=paste(type,marker,"fa", sep = ".") )
##
##
##marker="Inv_right"
##mask_info_inv_r = mask_info.all %>%
##  filter(V3 %in% c("inv"),
##         V2 > 6e6)
##inv_r_seq = inv2lt_seq
##colmask(inv_r_seq, invert=TRUE) <- mask_info_inv_r$adj_pos
##DNAStr_inv_r = as(inv_r_seq, "DNAStringSet")
##writeXStringSet(DNAStr_inv_r,file=paste(type,marker,"fa", sep = ".") )


