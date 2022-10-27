# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/jcbnunez/Shared_w_Alan/Local_enrrichment_df_summarize_pos_top.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)


## 1. Load object
load("~/Local_enrrichment_df_summarize_pos_top.Rdata")

#what is in this object?
names(Local_enrrichment_df_summarize_pos_top)



#category: Obs or Per, observed or permuted
#Symbol: Gene name
#Median_pos: Median position of the Gene in the GFF
#Analysis: All SNPs, all tranlated snos, Coding SNPs, Syn, NonSyn
#All_SNPs_tot: Number of total SNPs that fall in this category
#Median: This is the proportion of SNPs < P 5% / All SNPs. For the "Per" category this is the median over 100 perms. For "Obs" is the real value
#prop_l: output of prop_l = quantile(proportion, 0.01) for Perm
#prop_h: output of prop_l = quantile(proportion, 0.95) for Perm
#prop_sd: sd of permutations
#top_9X: Value of the top 95% of permutations this conects Obs and Per




#


gene <- as.data.table(Local_enrrichment_df_summarize_pos_top)
gene[,p:=1-pbinom(Median*All_SNPs_tot, All_SNPs_tot, .05)]

table(gene$p<.05, gene$category, gene$Analysis)

gene[category=="Obs"][Analysis=="NonSyn"][All_SNPs_tot>5][order(p)]
gene[Symbol=="Msp300"][All_SNPs_tot>5][order(p)]
gene[category=="Obs"][Analysis=="Translate"][All_SNPs_tot>=5][order(p)]
gene[Symbol=="CG17234"][All_SNPs_tot>5][order(p)]
gene[Symbol=="lectin-37Db"][All_SNPs_tot>5][order(p)]

gene[p==0, p:=1e-15]

p1 <-
ggplot(data=gene[All_SNPs_tot>5], aes(x=-log10(p), y=(All_SNPs_tot), color=Median, shape=as.factor(Median>top_9X))) +
geom_point() +
facet_grid(category~Analysis)

ggsave(p1, file="~/gene_enrichment.png")


unique(gene[All_SNPs_tot>5][p<1e-10][Median>top_9X[Analysis!="All"]]$Symbol)
unique(gene[Median>prop_h]$Symbol)


unique(gene[All_SNPs_tot>5][p<1e-10][Median>top_9X[Analysis!="All"]]$Symbol)



gene[Analysis=="Translate" & category=="Obs",pa:=p.adjust(p)]
gene[Analysis=="All" & category=="Obs",pa:=p.adjust(p)]

set <- gene[pa<.05][Median>top_9X][Analysis=="All"]$Symbol
universe <- gene[Analysis=="Translate"]$Symbol

write.csv(set, file="~/set.csv", quote=F, row.names=F)
write.csv(universe, file="~/universe.csv", quote=F, row.names=F)
