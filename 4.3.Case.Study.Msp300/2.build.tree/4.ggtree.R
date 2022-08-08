
library("treeio")
library("ggtree")
library("vroom")
library("TDbook")
library("adegenet")
library("FactoMineR")
library("tidyverse")


nwk <- "msp.tree/msp300.treefile"
tree <- read.tree(nwk)
#genotype <- as.data.frame(vroom("tree.metadat.txt"))
#genotype.tip = as.data.frame(genotype[,c("kar")])
#row.names(genotype.tip) = paste(genotype$tip, ".0", sep = "")
dna.seq <- fasta2DNAbin("Msp.300.haps.al.fasta", quiet=FALSE, chunkSize=10, snpOnly=FALSE)

bin.dna <- DNAbin2genind(dna.seq , pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100)

bin.dna@tab %>%
as.data.frame %>% 
.[,seq(from=1, to = dim(.)[2], by = 2 )] ->
simplified.bin.DNA

###
simplified.bin.DNA[-c(40:42),] %>% 
PCA(graph = F) ->
pca.obj

pca.obj$ind$coord %>%
as.data.frame %>% 
mutate(ind = rownames(.)) %>%
left_join(mutate(genotype, ind = paste(tip,".0", sep ="")  )) %>%  
ggplot(aes(
  x=Dim.2,
  y=Dim.3,
  color = kar,
  shape = pop
)) +
geom_point(size = 3) ->
pca.dna
ggsave(pca.dna, file = "pca.dna.pdf")

dimdesc(pca.obj) -> ddst
ddst$Dim.2 %>% as.data.frame %>% .[complete.cases(.),]

### pos 3844 in ALN

##
ggplot(tree, aes(x, y), 
branch.length="none"
) + 
geom_tree() + 
theme_tree() +
 geom_tiplab(size=2, align=TRUE, linesize=.5) + 
geom_nodelab(aes(x=branch, label=label), vjust=-.5, size=3) -> tree.plot


msaplot(tree.plot, dna.seq, 
window=c(998, 1008),  
offset = 5) -> dna.tree
ggsave(dna.tree, file = "dna.tree.tree.plot.pdf") 
