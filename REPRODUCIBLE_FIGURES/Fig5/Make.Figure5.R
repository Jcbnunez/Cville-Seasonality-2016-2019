## Figure 5
## 
library(tidyverse)
library(magrittr)
library(data.table)
library(gmodels)

### For Panel A
load("./out_count.Rdata")

### For pannel B
load("./co.enrich.files.Rdata")

## For panel C
dt.window =  get(load("./sliding.windowholm.Rdata"))

pheno <- readRDS("./wideform.fixed.phenotable.RDS")

### For panel D and E
pc_loading = readRDS("pca_pheno.data")
pcr.ag = readRDS("pca_loadings.data")
pca.info = load("pc.object.Rdata")

#### For Panels F/G
load(file = "deficiency.line.Rdata")
load(file = "phenotype.dt.Rdata")

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
              aes(x=group.general, y = N), alpha = 0.2, width = .2) +
  geom_point(shape = 23, 
             size = 5, 
             fill = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4") 
             ) +
  #facet_grid(cols = vars(inv), scales = "free")+
  theme_bw() + 
  ylim(0,20) +
  xlab("Phenotypic Categories") +
  ylab( "# of Sig. Phenotypes (< 0.05) ") 

  #theme(
  #  strip.text.x = element_text(size = 12),
  #  strip.text.y = element_text(size = 10),
  #  axis.title=element_text(size=12,face="bold"),
  #  axis.text.x = element_text( size = 12),
  #  axis.text.y = element_text(size = 12, face = "bold")) 

### Panel B

observedata %<>% filter(chr == "2L")
bounds %<>% filter(chr == "2L")

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
  facet_grid(~inv) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4"))  +
  scale_x_continuous(breaks = seq(0,3, by = 0.5)) +
  theme(legend.position = "none") -> panelB


#### Panel C
final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6,
         chr = "2L")

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
  geom_ribbon(data = dt.window[perm.st == "permutation"][chr == "2L"],
              aes(x=I(start/2+end/2)/1e6, ymax= upperbound, ymin = 0), fill="grey", linetype = "solid", alpha = 0.5) +
  geom_point(data=dt.window[perm.st == "observed"][chr == "2L"],
             aes(x=I(start/2+end/2)/1e6, y=avg.N, color = sig.v.per), #color="firebrick4",
             size = 2) +
  theme_bw()+
  theme(legend.position = "none") +
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos") -> PanelC

ggsave(PanelC, file =  "./new.model2.5.fdr.pdf", w = 7, h = 3)

### Panel D

#lim <- 1.8
p3 = ggplot(data=pc_loading) +
  #coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), size=0.1, arrow=arrow(length=unit(0.05,"cm"))) +
  geom_text(data=pc_loading[Dim.1<0], aes(x=x, y=y, label=pheno), size = 1.0, vjust=1, hjust=1) +
  geom_text(data=pc_loading[Dim.1>0], aes(x=x, y=y, label=pheno), size = 1.0, vjust=1, hjust=0) +
  geom_segment(aes(x=x, y=y, xend=Dim.1, yend=Dim.2), color="grey", alpha=0.5, size = 0.4) +
  xlim(-2.5,2.5) + ylim(-1.8,1.8) +
  theme(legend.position = "none") +
  theme_bw() + xlab("PCA Dimension 1 Loadings (8.7%)") + ylab("PCA Dimension 2 Loadings (7.8%)")

pc_loading %>%
  filter(Dim.1 < 0) %>%
  select(pheno)

pc_loading %>%
  filter(Dim.1 > 0) %>%
  select(pheno)

p1 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu1)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu1-1.96*se1, ymax=mu1+1.96*se1), width=.1) +
  geom_point() +
  theme_bw()  +
  xlab(NULL) +
  #coord_fixed() +
  ylab("PCA 1 Projection")

p2 <- ggplot(data=pcr.ag, aes(x=as.factor(gt.name), y=mu2)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu2-1.96*se2, ymax=mu2+1.96*se2), width=.1) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  #coord_fixed() +
  ylab("PCA 2 Projection")

p3+(p1/p2)  +
  plot_layout(guides = 'collect') -> panelD
####
###
###
###
###

ggplot(em.data, aes(x = f1.background, y = scaled.minute.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = inversion.st)) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  #facet_grid(.~f1.background, scales = "free_y")+
  theme_bw() + 
  #scale_color_manual(values = group.colors) +
  facet_grid("model.slope" ~ . ) +
  theme(legend.position = "none") +
  ylab("activity.slope") -> topDL
#part 2- startle duration

phenotype.dt[knockout == "k2a"] %>%
  group_by(inversion.st, f1.background,  knockout,# inversion.st 
  ) %>%
  summarise(mean = ci(sr.length)[1],
            uci = ci(sr.length)[2],
            lci = ci(sr.length)[3]
  ) %>%
  ggplot(aes(
    x=f1.background,
    y=mean,
    ymin=lci,
    ymax=uci,
    color =  inversion.st
  )) +
  # scale_color_manual(values = group.colors) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  facet_grid(knockout~., scales = "free_y")+
  xlab("gene.deficiecy.status") +
  ylab("startle.length") +
  theme(legend.position = "none") +
  theme_bw() -> bottomDL

(topDL/bottomDL) -> panelE

(panelA+panelB+PanelC)/(panelD+panelE)  +
  plot_layout(guides = 'collect')

ggsave(panelA, file = "panelA.pdf", w = 1.5, h = 1.5)
ggsave(panelB, file = "panelB.pdf", w = 3.0, h = 1.5)
ggsave(PanelC, file = "panelC.pdf", w = 4.0, h = 1.5)
ggsave(panelD, file = "panelD.pdf", w = 7.0, h = 2.5)
ggsave(panelE, file = "panelE.pdf", w = 5.0, h = 2.5)


