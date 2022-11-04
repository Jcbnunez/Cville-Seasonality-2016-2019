#system("cp /project/berglandlab/Yang_Adam/reference_files/dgrp.gds .")
wd =  paste0("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/13.Phenotype_Analysis/DATA")
#wd = "/scratch/bal7cg/Inverson_project/"
setwd( wd )
### libraries
library(data.table)
library(ggplot2)
library(SeqArray)
library("Hmisc")
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(patchwork)

### load
pheno <- readRDS("./wideform.fixed.phenotable.RDS")
dim(pheno)

pl <- melt(pheno, id.vars="ral_id")
pl.ag <- pl[,list(missing=sum(is.na(value)),
                  mean=mean(value, na.rm=T),
                  sd=sd(value, na.rm=T)), list(variable)]


pl <- merge(pl, pl.ag, by="variable")
pl[,x:=value]
pl[is.na(value), x:=mean]
pl[,x_norm:=(x-mean)/sd]
pl[,ral_id:=paste("line_", ral_id, sep="")]

#load in new phenos
newphenos = readRDS("phenossaved")
newphenos = data.table(
  pheno = newphenos
)
colnames(newphenos)[1] = "pheno"
newphenos[,pheno:= tstrsplit(pheno, "\\.")[[1]]]

### inversion
inv <- fread("./inversion.csv")
names(inv)[1] = "ral_id"
inv$ral_id = gsub("DGRP", "line", inv$ral_id)

### principal components
#foreach(unique(pl$pos))%do%{

#pw <-dcast(pl[pos==5192177][variable%in%target$pheno], ral_id+gt~variable, value.var="x_norm")
pw = dcast(pl[variable %in% newphenos$pheno],ral_id ~ variable, value.var = "x_norm")
#w <-dcast(pl, ral_id+gt~variable, value.var="x_norm")

pc <- PCA(pw[,-c("ral_id"), with=F], scale.unit=T, graph=F)
#pc = PCA(pheno, scale.unit = T, graph = F)
plot.PCA(pc)
#viz <- fviz_pca_var(pc, select.var = list(name = NULL, cos2 = NULL, contrib = NULL), repel=F) + theme(text=element_text(size=1))


pc_loading <- as.data.table(pc$var$coord)
pc_loading[,pheno:=rownames(pc$var$coord)]
step=.1
### top left quadrant
pc_loading[Dim.1<0 & Dim.2>0, x:=-.4]
pc_loading[Dim.1<0 & Dim.2>0, 
           y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2>0])[1], by=step)]

### bottom left quadrant
pc_loading[Dim.1<0 & Dim.2<0, x:=-.25]
pc_loading[Dim.1<0 & Dim.2<0, 
           y:=-1 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2<0])[1], by=step)]

### bottom right quadrant
pc_loading[Dim.1>0 & Dim.2<0, x:=.25]
pc_loading[Dim.1>0 & Dim.2<0, 
           y:=-.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.2<0])[1], by=step)]

### top right quadrant lower on Dim1 (TOP CLUSTER)
pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, 
           x:=.05]
pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, 
           y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0])[1], by=step)]

### top right quadrant upper on Dim1
pc_loading[Dim.1>.25 & Dim.2>0, x:=.8]
pc_loading[Dim.1>.25 & Dim.2>0, 
           y:=.65 - seq(from=0, length.out=dim(pc_loading[Dim.1>.25 & Dim.2>0])[1], by=step)]

### indiviuals PCA data table
pcr <- cbind(pw[,c("ral_id"), with=F], as.data.table(pc$ind$coord))
pcr <- merge(pcr, inv[`In(2L)t` != "INV/ST"], by="ral_id") #remove heterozygotes for now

# summary(lm(Dim.1~gt, pcr))
# summary(lm(Dim.2~gt, pcr))
# summary(lm(Dim.3~gt, pcr))

summary(lm(Dim.1~`In(2L)t`, pcr))
summary(lm(Dim.2~`In(2L)t`, pcr))
summary(lm(Dim.3~`In(2L)t`, pcr))




pcr.ag <- pcr[!is.na(`In(2L)t`),list(mu1=mean(Dim.1), se1=sd(Dim.1)/sqrt(length(Dim.1)),
                              mu2=mean(Dim.2), se2=sd(Dim.2)/sqrt(length(Dim.2)),
                              mu3=mean(Dim.3), se3=sd(Dim.3)/sqrt(length(Dim.3))),
              list(`In(2L)t`)]
pcr.ag[`In(2L)t`== "INV", gt.name:="Inverted"]
pcr.ag[`In(2L)t`== "ST", gt.name:="Standard"]

### loadings plot
pcl <- pc$var$coord
pclw <- as.data.table(melt(pcl))
pclw[,dim:=as.numeric(tstrsplit(Var2, "\\.")[[2]])]
#ggplot(data=pclw, aes(x=dim, y=value, group=Var1)) + geom_line() +
#geom_text(data=pclw[dim==1], aes(x=1, y=value, label=Var1), size=2, hjust="right") +
#xlim(-1, 5)
saveRDS(pc_loading, "pca_pheno.data")
saveRDS(pcr.ag, "pca_loadings.data")
p1 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu1)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu1-1.96*se1, ymax=mu1+1.96*se1), width=.1) +
  geom_point() +
  theme_bw()  +
  xlab(NULL) +
  ylab("PCA Dimension 1\nProjection")

ggsave(p1, file = "p1.pdf")

p2 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu2)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu2-1.96*se2, ymax=mu2+1.96*se2), width=.1) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  ylab("PCA Dimension 2\nProjection")


lim <- 1.8
 p3 = ggplot(data=pc_loading) +
  coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data=pc_loading[Dim.1<0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=1) +
  geom_text(data=pc_loading[Dim.1>0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=0) +
  
  geom_segment(aes(x=x, y=y, xend=Dim.1, yend=Dim.2), color="grey", alpha=0.5) +
  xlim(-lim,lim) + ylim(-lim,lim) +
  theme_bw() + xlab("PCA Dimension 1 Loadings (12.5%)") + ylab("PCA Dimension 2 Loadings (7.8%)")


