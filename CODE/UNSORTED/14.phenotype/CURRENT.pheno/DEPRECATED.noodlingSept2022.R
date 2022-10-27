system("cp /project/berglandlab/Yang_Adam/reference_files/dgrp.gds .")

### libraries
  library(data.table)
  library(ggplot2)
  library(SeqArray)
  library("Hmisc")
  library(FactoMineR)
  library(factoextra)

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

### load GLM
  glm <- load("/Users/alanbergland/glm_flight/glm.out.VA_ch_0.Rdata")
  best <- glm.out[mod=="aveTemp+year_factor"][chr=="2L"][pos>=(5192177-10000)][pos<=(5192177+10000)]
  #best <- glm.out[mod=="aveTemp+year_factor"][chr=="2L"][pos==5192177]
  best <- best[!is.na(rnp.clean)]
  setkey(snp.dt, chr, pos)
  setkey(best, chr, pos)
  best <- snp.dt[J(best[!is.na(rnp.clean)])]

### open DGRP genotype data
  genofile <- seqOpen("~/dgrp.gds")
  snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                       chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"))

###
  ### TSP: 5192177
  ### best in region: 5191923
  
sigo <- foreach(v.i=best$variant.id, .combine="rbind")%do%{
  #pos.i <- 5192177
  # v.i <- 247978

  seqSetFilter(genofile, variant.id=v.i)

  genos <- data.table(ral_id=seqGetData(genofile, "sample.id"),
                      variant.id=rep(seqGetData(genofile, "variant.id"), each=205),
                      pos=rep(seqGetData(genofile, "position"), each=205),
                      gt=expand.grid(seqGetData(genofile, "$dosage")[,1]))
  setnames(genos, "gt.Var1", "gt")

  ### merge
    pl.use <- merge(pl, genos, by="ral_id", allow.cartesian=T)

  ### principal components
    pw <-dcast(pl.use[variable%in%target$pheno], ral_id+gt~variable, value.var="x_norm")
    #w <-dcast(pl, ral_id+gt~variable, value.var="x_norm")

    pc <- PCA(pw[,-c("ral_id", "gt"), with=F], scale.unit=T, graph=F)
    #plot.PCA(pc)
    #viz <- fviz_pca_var(pc, select.var = list(name = NULL, cos2 = NULL, contrib = NULL), repel=F) + theme(text=element_text(size=1))

    ### indiviuals PCA data table
      pcr <- cbind(pw[,c("ral_id", "gt"), with=F], as.data.table(pc$ind$coord))
      pcr <- merge(pcr, inv, by="ral_id")

    t1 <- (lm(Dim.1~gt, pcr[gt!=1]))
    t2 <- (lm(Dim.2~gt, pcr[gt!=1]))
    t3 <- (lm(Dim.3~gt, pcr[gt!=1]))

    data.table(pos=best[variant.id==v.i]$pos, pc=c(1,2,3), rnp.clean=best[variant.id==v.i]$rnp.clean,
                maf=sum(pcr[gt!=1]$gt/2)/length(pcr[gt!=1]$gt),
                p=c(summary(t1)$coef[2,4], summary(t2)$coef[2,4], summary(t3)$coef[2,4]))
  }
  fisher.test(table(sigo$p<.0005, sigo$rnp.clean<.05))
  table(sigo$p<.005, sigo$rnp.clean<.05, sigo$pc)
  sigo[pc==1][which.min(p)]
  sigo[pc==1][which.min(rnp.clean)]
  sigo[pos==5192177]

    pc_loading <- as.data.table(pc$var$coord)
    pc_loading[,pheno:=rownames(pc$var$coord)]
    step=.1
    ### top left quadrant
      pc_loading[Dim.1<0 & Dim.2>0, x:=-.4]
      pc_loading[Dim.1<0 & Dim.2>0, y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2>0])[1], by=step)]

    ### bottom left quadrant
      pc_loading[Dim.1<0 & Dim.2<0, x:=-.25]
      pc_loading[Dim.1<0 & Dim.2<0, y:=-1 - seq(from=0, length.out=dim(pc_loading[Dim.1<0 & Dim.2<0])[1], by=step)]

    ### bottom right quadrant
      pc_loading[Dim.1>0 & Dim.2<0, x:=.25]
      pc_loading[Dim.1>0 & Dim.2<0, y:=-.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.2<0])[1], by=step)]

    ### top right quadrant lower on Dim1 (TOP CLUSTER)
      pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, x:=.05]
      pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0, y:=1.5 - seq(from=0, length.out=dim(pc_loading[Dim.1>0 & Dim.1 <=0.25 & Dim.2>0])[1], by=step)]

    ### top right quadrant upper on Dim1
      pc_loading[Dim.1>.25 & Dim.2>0, x:=.8]
      pc_loading[Dim.1>.25 & Dim.2>0, y:=.65 - seq(from=0, length.out=dim(pc_loading[Dim.1>.25 & Dim.2>0])[1], by=step)]



  pcr.ag <- pcr[!is.na(gt),list(mu1=mean(Dim.1), se1=sd(Dim.1)/sqrt(length(Dim.1)),
                                mu2=mean(Dim.2), se2=sd(Dim.2)/sqrt(length(Dim.2)),
                                mu3=mean(Dim.3), se3=sd(Dim.3)/sqrt(length(Dim.3))),
                          list(gt)]
  pcr.ag[gt==0, gt.name:="Cold allele"]
  pcr.ag[gt==2, gt.name:="Hot allele"]

### loadings plot
  pcl <- pc$var$coord
  pclw <- as.data.table(melt(pcl))
  pclw[,dim:=as.numeric(tstrsplit(Var2, "\\.")[[2]])]
  #ggplot(data=pclw, aes(x=dim, y=value, group=Var1)) + geom_line() +
  #geom_text(data=pclw[dim==1], aes(x=1, y=value, label=Var1), size=2, hjust="right") +
  #xlim(-1, 5)

  p1 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu1)) +
        geom_errorbar(aes(x=as.factor(gt.name), ymin=mu1-1.96*se1, ymax=mu1+1.96*se1), width=.1) +
        geom_point() +
        theme_bw()  +
        xlab(NULL) +
        ylab("PCA Dimension 1\nProjection")

  p2 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu2)) +
        geom_errorbar(aes(x=as.factor(gt.name), ymin=mu2-1.96*se2, ymax=mu2+1.96*se2), width=.1) +
        geom_point() +
        theme_bw() +
        xlab(NULL) +
        ylab("PCA Dimension 2\nProjection")


  lim <- 1.8
  viz <- ggplot(data=pc_loading) +
         coord_equal() +
         geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm"))) +
         geom_text(data=pc_loading[Dim.1<0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=1) +
         geom_text(data=pc_loading[Dim.1>0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=0) +

         geom_segment(aes(x=x, y=y, xend=Dim.1, yend=Dim.2), color="grey", alpha=0.5) +
         xlim(-lim,lim) + ylim(-lim,lim) +
         theme_bw() + xlab("PCA Dimension 1 Loadings (12.5%)") + ylab("PCA Dimension 2 Loadings (7.8%)")



  layout <- "
  AA#
  AAB
  AAC
  AA#"

  mega <- viz + p1 + p2 +  plot_layout(design=layout)

  ggsave(mega, file="~/mega3.pdf", h=8, w=8)







### cmds
cmd <- cmdscale(dist(pw[,-c("ral_id", "gt"), with=F]))
cmd <- cbind(as.data.table(cmd), pw[,"ral_id", with=F])
cmd <- merge(cmd, inv, by="ral_id")
cmd <- cbind(pw[,c("ral_id", "gt"), with=F], cmd)
summary(lm(V1~gt, cmd))
summary(lm(Dim.2~gt, pcr))
summary(lm(Dim.3~gt, pcr))



### tsne
library(uwot)
umap_cereals_num <- umap(t(pw[,-c("ral_id", "gt"), with=F]),
  n_neighbors = 15,
  min_dist = 1, spread = 5
)

umap_cereals_num <- data.table(
  UMAP1 = umap_cereals_num[, 1],
  UMAP2 = umap_cereals_num[, 2],
  label = row.names(umap_cereals_num)
)

ggplot(umap_cereals_num, aes(
  x = UMAP1, y = UMAP2,
  label = label
)) +
  geom_point() +
  ggrepel::geom_text_repel(cex = 2.5)





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
