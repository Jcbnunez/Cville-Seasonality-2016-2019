system("scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/Yang_Adam/reference_files/dgrp.gds ~/.")

### libraries
  library(data.table)
  library(ggplot2)
  library(SeqArray)
  library("Hmisc")


### load
  pheno <- readRDS("~/wideform.phenotypedata.RDS")
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

### target phenotypes
  target <- fread("/Users/alanbergland/hotandcold.phenos.txt")
  target <- data.table(pheno=c(target$hot, target$cold), polar=rep(c("hot", "cold"), each=dim(target)[1]))
  target <- na.omit(target)

### inversion
  inv <- fread("~/inversion.csv")

### open DGRP genotype data
  genofile <- seqOpen("~/dgrp.gds")
  snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                       chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"))

  seqSetFilter(genofile, variant.id=snp.dt[pos==5192177]$variant.id)
  genos <- data.table(ral_id=seqGetData(genofile, "sample.id"),
                      gt=seqGetData(genofile, "$dosage"))
  setnames(genos, "gt.V1", "gt")

### merge
  pl <- merge(pl, genos, by="ral_id")

### principal components
  pw <-dcast(pl[variable%in%target$pheno], ral_id+gt~variable, value.var="x_norm")
  #w <-dcast(pl, ral_id+gt~variable, value.var="x_norm")

  pc <- PCA(pw[,-c("ral_id", "gt"), with=F])
  plot.PCA(pc)
  viz <- fviz_pca_var(pc, select.var = list(name = NULL, cos2 = NULL, contrib = 15), repel=F)

  ggplot(data=pc$ind$coord, aes(x=Dim.1, y=Dim.2)) + geom_point()

  pcr <- cbind(pw[,c("ral_id", "gt"), with=F], as.data.table(pc$ind$coord))
  pcr <- merge(pcr, inv, by="ral_id")

  summary(lm(Dim.1~gt, pcr))
  summary(lm(Dim.2~gt, pcr))
  summary(lm(Dim.3~gt, pcr))


  pcr.ag <- pcr[!is.na(gt),list(mu1=mean(Dim.1), se1=sd(Dim.1)/sqrt(length(Dim.1)),
                                mu2=mean(Dim.2), se2=sd(Dim.2)/sqrt(length(Dim.2)),
                                mu3=mean(Dim.3), se3=sd(Dim.3)/sqrt(length(Dim.3))),
                          list(gt)]
  pcr.ag[gt==0, gt.name:="Cold allele"]
  pcr.ag[gt==2, gt.name:="Hot allele"]

  p1 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu1)) +
        geom_errorbar(aes(x=as.factor(gt.name), ymin=mu1-1.96*se1, ymax=mu1+1.96*se1), width=.1) +
        geom_point() +
        theme_bw() + coord_flip()
  p2 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu2)) +
        geom_errorbar(aes(x=as.factor(gt.name), ymin=mu2-1.96*se2, ymax=mu2+1.96*se2), width=.1) +
        geom_point() +
        theme_bw() + coord_flip()
  p3 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu3)) +
        geom_errorbar(aes(x=as.factor(gt.name), ymin=mu3-1.96*se3, ymax=mu3+1.96*se3), width=.1) +
        geom_point() +
        theme_bw() + coord_flip()

  layout <- "
  AB
  CD"

  mega <- viz + p1 + p2 + p3 + plot_layout(design=layout)

  ggsave(mega, file="~/mega2.pdf", h=8, w=8)


### some basic correlation matrices
  pw <-dcast(pl[variable%in%target$pheno][gt==0], ral_id~variable, value.var="x_norm")
  pwc <- rcorr(as.matrix(pw[,-1]))
  pwc.dt <- cbind(data.table(pheno=names(pw)[-1]), as.data.table(pwc$r))
  heatmap(pwc$r)
  diag(pwc$P) <- 1

  p0 <-
  ggcorrplot(corr=(pwc$r), p.mat=pwc$P,
              hc.order = T, outline.col = "white") +
  theme(axis.text.y =  element_text(size=4),
        axis.title.y = element_text(size=4),
        axis.ticks.y = element_blank(),
         axis.text.x = element_text(size=4),
        axis.title.x = element_text(size=4),
        axis.ticks.x = element_blank(),
        legend.position = "top")



  pw <-dcast(pl[variable%in%target$pheno][gt==2], ral_id~variable, value.var="x_norm")
  pwc <- rcorr(as.matrix(pw[,-1]))
  pwc.dt <- cbind(data.table(pheno=names(pw)[-1]), as.data.table(pwc$r))
  heatmap(pwc$r)
  diag(pwc$P) <- 1

  p1 <-
  ggcorrplot(corr=(pwc$r), p.mat=pwc$P,
              hc.order = T, outline.col = "white") +
  theme(axis.text.y =  element_text(size=4),
        axis.title.y = element_text(size=4),
        axis.ticks.y = element_blank(),
         axis.text.x = element_text(size=4),
        axis.title.x = element_text(size=4),
        axis.ticks.x = element_blank(),
        legend.position = "top")

  p0 + p1






















###### Defunc
###
  pw <-dcast(pl[variable%in%target$pheno], ral_id+gt~variable, value.var="x_norm")
  #w <-dcast(pl, ral_id+gt~variable, value.var="x_norm")


  pc <- PCA((pw[,-c("ral_id", "gt"), with=F]))
ddpc <- dimdesc(pc)

tmp1 <- rbind(data.table(pheno=rownames(ddpc$Dim.1$quanti), ddpc$Dim.1$quanti))
tmp1[,x:=c(1:dim(tmp1)[1])]

tmp2 <- rbind(data.table(pheno=rownames(ddpc$Dim.2$quanti), ddpc$Dim.2$quanti))
tmp2[,x:=c(1:dim(tmp2)[1])]

tmp <- merge(tmp1, tmp2, by="pheno")
tmp[,pheno:=factor(pheno, levels=pheno[x.x])]

ggplot(data=tmp) +
geom_point(aes(x=correlation.x, y=x.x), color="red") +
geom_point(aes(x=correlation.y, y=x.x), color="blue")

nb.dim <- estim_ncp(pw[,-c("ral_id", "gt")],scale=TRUE)


res.mca <- MCA(as.data.table(), quanti.sup=1)

  mds <- cmdscale(dist(t(pw[,-c("ral_id", "gt"), with=F])), k=2)
  pcr <- cbind(pw[,c("ral_id", "gt"), with=F], as.data.table(mds))
  pcr <- merge(pcr, inv, by="ral_id")








    pc <- prcomp((pw[,-c("ral_id", "gt"), with=F]))
    pw.pc <- cbind(pw[,c("ral_id", "gt"), with=F],
                    pc$x[,1:10])
    pcr <- as.data.table(pc$rotation)
    pcr[,pheno:=rownames(pc$rotation)]
    pcr[order(-abs(PC1))][,c("pheno", "PC1", "PC2"), with=F][1:10]

    summary(lm(PC1~gt, pw.pc))
    summary(lm(PC2~gt, pw.pc))
    summary(lm(PC3~gt, pw.pc))
    summary(lm(PC4~gt, pw.pc))

  pw <-dcast(pl, ral_id~variable, value.var="x_norm")


library("Hmisc")
  pwc <- rcorr(as.matrix(pw[,-1]))
  pwc.dt <- cbind(data.table(pheno=names(pw)[-1]), as.data.table(pwc$r))

  # Cluster & heatmap on otter data
  # Jeff Oliver
  # jcoliver@email.arizona.edu
  # 2017-08-15

  ################################################################################

  # Read in data

  # Scale each measurement (independently) to have a mean of 0 and variance of 1
  otter_scaled <- as.data.frame(pwc.dt)
  #otter_scaled[, -1] <- scale(otter_scaled[, -1])

  # Run clustering
  otter_matrix <- as.matrix(otter_scaled[,-1])
  rownames(otter_matrix) <- otter_scaled$pheno
  otter_dendro <- as.dendrogram(hclust(d = dist(x = otter_matrix)))

  # Create dendrogram plot
  dendro_plot <- ggdendrogram(data = otter_dendro, rotate = TRUE) +
    theme(axis.text.y = element_text(size = 6))

  # Heatmap

  # Data wrangling
  otter_long <- pivot_longer(data = otter_scaled,
                             cols = -c(pheno),
                             names_to = "measurement",
                             values_to = "value")
  # Extract the order of the tips in the dendrogram
  otter_order <- order.dendrogram(otter_dendro)
  # Order the levels according to their position in the cluster
  otter_long$accession <- factor(x = otter_long$pheno,
                                 levels = otter_scaled$pheno[otter_order],
                                 ordered = TRUE)

  # Create heatmap plot
  heatmap_plot <- ggplot(data = otter_long, aes(x = measurement, y = pheno)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top")

  # All together

  p1 <- grid.newpage()
  print(heatmap_plot,
        vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
  print(dendro_plot,
        vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.92))

ggsave(p1, file="~/p1.pdf", h=8, w=8)



heatmap(
  as.matrix(pwc$r), Rowv=NA,
  Colv=as.dendrogram(hclust(dist(t(as.matrix(pwc$r)))))
)



pw <-dcast(pl, ral_id~variable, value.var="x_norm")
p.mat <- cor_pmat(pw[,-1])

p1 <- ggcorrplot(p.mat, hc.order = F, outline.col = "white") +
theme(axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
       axis.text.y=element_blank(),
       axis.ticks.y=element_blank()
       )


       p2 <- ggcorrplot(p.mat, hc.order = T, outline.col = "white") +
       theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()
              )


library(patchwork)
p1 + p2



library(plotly)

rownames(pwc$r) <- NULL
colnames(pwc$r) <- NULL
p1 <- heatmaply_cor(
  (pwc$r))



  +
  theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()
         )

ggsave(p1, file="~/test.pdf")


  pwd <- dist(pw)
  pwh <- hclust(pwd)
  pwc <- cor(pw)
  ggplot(data=pl, aes(x=ral_id, y=variable, fill=value)) + geom_tile()


  heatmap(as.matrix(pw)[,-1])
