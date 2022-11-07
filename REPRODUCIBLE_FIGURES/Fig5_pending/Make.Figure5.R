## Figure 5
## 


load("./out_count.Rdata")
load("./co.enrich.files.Rdata")
dt.window =  get(load("./sliding.windowholm.Rdata"))
pheno <- readRDS("./wideform.fixed.phenotable.RDS")
pc_loading = readRDS("pca_pheno.data")
pcr.ag = readRDS("pca_loadings.data")

### Panel A
#shorten names
out_count$group.general = gsub("behavior", "B.", out_count$group.general)
out_count$group.general = gsub("stress-resistance", "SR.", out_count$group.general)
out_count$group.general = gsub("life-history", "LH.", out_count$group.general)
out_count$group.general = gsub("morphology", "M.", out_count$group.general)
#add inversion status to chr
out_count$inv = out_count$chr
out_count$inv = gsub("2L", "2L Inversion Models", out_count$inv)
out_count$inv = gsub("2R", "2R Inversion Models", out_count$inv)
out_count$inv = gsub("3L", "3L Inversion Models", out_count$inv)
out_count$inv = gsub("3R", "3R Inversion Models", out_count$inv)


#try graph again
panelA = ggplot(data = out_count[perm.status == "observed"], aes(x = group.general, y = N)) +
  geom_boxplot(width = 0.5, data = out_count[perm.status == "permutation"], 
               aes(x=group.general, y = N), outlier.colour = "NA") +
  geom_jitter(data = out_count[perm.status == "permutation"], 
              aes(x=group.general, y = N), alpha = 0.08, width = .2) +
  geom_point(shape = 23, 
             size = 5, 
             fill = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4") 
             ) +
  #facet_grid(cols = vars(inv), scales = "free")+
  theme_bw() + 
  xlab("Phenotypic Categories") +
  ylab( "# of Sig. Phenotypes (< 0.05) ") #+
  #theme(
  #  strip.text.x = element_text(size = 12),
  #  strip.text.y = element_text(size = 10),
  #  axis.title=element_text(size=12,face="bold"),
  #  axis.text.x = element_text( size = 12),
  #  axis.text.y = element_text(size = 12, face = "bold")) 

### Panel B

bounds$inv = gsub(FALSE, "Not-Inverted", bounds$inv)
bounds$inv = gsub(TRUE, "Inverted", bounds$inv)
observedata$inv = gsub(FALSE, "Not-Inverted", observedata$inv)
observedata$inv = gsub(TRUE, "Inverted", observedata$inv)

ggplot() +
  geom_point(data = observedata, aes(x = or, y = prop, color = group.general)) +
  geom_abline(data = bounds, aes(intercept = mean.lower.prop, slope = 0), color = "black", linetype = "dashed")+
  geom_abline(data = bounds, aes(intercept = mean.upper.prop, slope = 0), color = "black", linetype = "dashed")+
  geom_vline(data = bounds, aes(xintercept = mean.lower.or), color = "black", linetype = "dashed")+
  geom_vline(data = bounds, aes(xintercept = mean.upper.or), color = "black", linetype = "dashed")+
  theme_bw()+
  xlab("coenrichment between GWAS/GLM") +
  ylab(" codirectionality in GWAS/GLM") +
  facet_grid(cols = vars(chr),rows = vars(inv)) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4"))  +
  scale_x_continuous(breaks = seq(0,3, by = 0.5))

+
  theme(
    #legend.position = c(0.875, 2.05),
    legend.position = "none",
    legend.key.size = unit(.4, "cm") , 
    legend.text = element_text(size = 8),
    axis.text.x = element_text( size = 7, face = "bold"),
  ) 



#### Panel C
ggplot() +
  #geom_rect(data=outlier_haplowins, aes(xmin=start/1e6, xmax=end/1e6, ymin=-1, ymax=100),fill ="lightgoldenrod1", alpha=.5) +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = 0, ymax = 40), 
            alpha = 0.7, fill = "gold") +
  #geom_line(data = dt.window[perm.st == "permutation"],
  #           aes(x=I(start/2+end/2)/1e6, y=avg.N), color="black") +
  #geom_line(data = dt.window[perm.st == "permutation"],
  #           aes(x=I(start/2+end/2)/1e6, y=lowerbound), color="black", linetype = "dashed") +
  geom_vline(xintercept = 2225744/1e6) +
  geom_vline(xintercept = 13154180/1e6) +
  geom_ribbon(data = dt.window[perm.st == "permutation"],
              aes(x=I(start/2+end/2)/1e6, ymax= upperbound, ymin = 0), fill="grey", linetype = "solid", alpha = 0.5) +
  geom_point(data=dt.window[perm.st == "observed"],
             aes(x=I(start/2+end/2)/1e6, y=avg.N, color = sig.v.per), #color="firebrick4",
             size = 2) +
  
  theme_bw()+
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos") -> g1

ggsave(g1, file =  "./new.model2.5.fdr.pdf", w = 7, h = 3)




### Panel D

lim <- 1.8
p3 = ggplot(data=pc_loading) +
  coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data=pc_loading[Dim.1<0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=1) +
  geom_text(data=pc_loading[Dim.1>0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=0) +
  geom_segment(aes(x=x, y=y, xend=Dim.1, yend=Dim.2), color="grey", alpha=0.5) +
  xlim(-lim,lim) + ylim(-lim,lim) +
  theme(legend.position = "none") +
  theme_bw() + xlab("PCA Dimension 1 Loadings (17.3%)") + ylab("PCA Dimension 2 Loadings (13.2%)")

p1 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu1)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu1-1.96*se1, ymax=mu1+1.96*se1), width=.1) +
  geom_point() +
  theme_bw()  +
  xlab(NULL) +
  coord_fixed() +
  ylab("PCA 1 Projection")

p2 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu2)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu2-1.96*se2, ymax=mu2+1.96*se2), width=.1) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  coord_fixed() +
  ylab("PCA 2 Projection")

p3+(p1/p2)
####

