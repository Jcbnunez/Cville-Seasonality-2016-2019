library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(adegenet)
library(data.table)
library(reshape2)
library(pegas)
library(vcfR)

### Prepare filtering parameters
### Finding maf
#vcf_maf = maf(vcf, element = 2)
#vcf_maf %>%
#  as.data.frame() %>%
#  mutate(row_id = 1:dim(.)[1]) ->
#  maf_master 
#
#maf_master%>%
#  filter(Frequency > 0.05) %>%
#  .$row_id -> maf_filter
#
## Make DNAbin




my_dnabin1 <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE)
my_dnabin1

glm_ld_markers_gen = DNAbin2genind(x = my_dnabin1, polyThres = 0.05)
glm_ld_markers_gen@tab %>%
  as.data.frame -> snp_dat_table
snp_dat_table[,seq(from=1, to=dim(snp_dat_table)[2], by=2)] -> snp_dat_table_loc
dim(snp_dat_table_loc)

set.seed(123)
res <- hcut(snp_dat_table_loc, k = 4, stand = TRUE)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF"
                       # ,"#E7B800", "#FC4E07"
          )) -> cluster_tree
ggsave(cluster_tree, file = "cluster_tree.pdf")

######

##finding invs only
for(karyo in 1:4  ){
  
  pdf(paste(karyo,"MSA.inv.inv2lt.pdf", sep = "."))
  ape::image.DNAbin(my_dnabin1[names(res$cluster[res$cluster == karyo])
                               ,maf_filter], show.labels = F)
  dev.off()
}


snp_dat_table %>%
  PCA(., 
      graph = FALSE,
      scale.unit = F,
      ncp=5 ) -> pca_fig_std

pca_fig_std$ind$coord %>% 
  as.data.frame %>% 
  mutate(samp_hap =  rownames(.)) %>% 
  mutate(samp = gsub("_0$|_1$", "", samp_hap)  ) %>% 
  left_join(., metadata_tab) %>%  
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = inv_st
  )) +
  geom_point() -> pca_gg_figure

ggsave(pca_gg_figure, file = "pca_gg_figure.pdf")









dataHaplo<-haplotype(dna_seq)
dataHaplo

d <- dist.dna(dataHaplo, "N")
nt <- rmst(d, quiet = TRUE)

pdf("raw_hap.pdf")
setHaploNetOptions(haplotype.inner.color = "#CCCC4D",
                   haplotype.outer.color = "#CCCC4D",
                   show.mutation = 3, labels = FALSE)
plot(nt, size = 2)
dev.off()





#######
#######
#######
#######
#######
#######
#######
#######extra


glm_ld_markers <- fasta2DNAbin("./Homozyg_only.GlmInv.fa", 
                               quiet=FALSE, chunkSize=10, snpOnly=FALSE)
glm_ld_markers_gen = DNAbin2genind(x = glm_ld_markers, polyThres = 0.00)
glm_ld_markers_gen@tab %>%
  as.data.frame -> snp_dat_table

snp_dat_table %>% dim

##glm_ld_markers.h <- fasta2DNAbin("/project/berglandlab/jcbnunez/Cville_2lt_haps/Het_only_OutlierInvs_rm.fa", 
##                               quiet=FALSE, chunkSize=10, snpOnly=FALSE)
##glm_ld_markers.h_gen = DNAbin2genind(x = glm_ld_markers.h, polyThres = 0.00)
##glm_ld_markers.h_gen@tab %>%
##  as.data.frame -> het_dat_table
### join
#rbind(homozyg_dat_table, het_dat_table) -> snp_dat_table

#Do an additional filtering step (this is related to the way that R load the data) 
snp_dat_table[,seq(from=1, to=dim(snp_dat_table)[2], by=2)] -> snp_dat_table_loc
dim(snp_dat_table_loc)

snp_dat_table_loc %>%
  as.data.frame %>%
  mutate(samp=rownames(.)) %>% 
  separate(samp, into = c("samp","haplo"), sep = "\\.") %>%
  left_join(metadata_tab)  %>%
  filter(inv_st %in% c("Std", "Inv"))  ->
  snp_dat_table_loc_meta

snp_dat_table_loc_meta %>%
  dplyr::select(-samp, -haplo, -inv_st) %>% 
  PCA(., 
      graph = FALSE,
      scale.unit = F,
      ncp=5 ) -> pca_fig_std

pca_fig_std$ind$coord %>%
  as.data.frame %>%
  cbind(., snp_dat_table_loc_meta[,c("samp", "haplo", "inv_st")]) %>% 
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color = haplo,
    shape = inv_st
  )) +
  geom_point() -> pca_gg_figure

ggsave(pca_gg_figure, file = "pca_gg_figure.pdf")


#### contriburions
## import snps
pos_adj_scalar = 2051600

mask_info.all <- fread("./Retain_loci_metadat.txt")
names(mask_info.all) = c("chr", "pos", "type")
mask_info.all %<>%
  mutate(adj_pos = pos-pos_adj_scalar+1)

mask_info.all %>%
  .$type %>% table


pca_fig_std$var$contrib %>%
  as.data.frame() %>% 
  cbind(mask_info.all) %>% 
  melt(id = c("chr", "pos", "type","adj_pos")) %>% 
  ggplot(aes(x=pos,
             y=value,
  )) +
  geom_point() +
  facet_wrap(~variable)->
  contrib_plot

ggsave(contrib_plot, file = "contrib_plot.pdf")
