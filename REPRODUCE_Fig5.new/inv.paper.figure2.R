###Inversion paper figure
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(corrplot)
library(patchwork)
library(ggrepel)
library(reshape2)
library(gmodels)
#install.packages('doParallel', repos='http://cran.us.r-project.org')
#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(10)
#setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
#load in data
load("figuredataOctober.Rdata")
# setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/May_2022_objects/")
# #load in data 
# #broad coenrichment data
# load("co.enrich.files")
# setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/Aug_2022_objects/")
# #karyotype modeling data
# 
# out_count = readRDS("karyotype.modeling.RDS")
# 
# 
# #slding window data
# # groupby.perm = readRDS("newmodelplotdata.RDS")
# dt.window = readRDS("sliding.windowholm.RDS")
# 
# # phenotype pca data
# pc_loading = readRDS("pca_pheno.data")
# #pca loading- inversion data
# pcr.ag = readRDS("pca_loadings.data")
# #deficiency data: emm counts
# em.data = readRDS("trends.startle.decay")
# 
# #deficiency data: startle length
# phenotype.dt = readRDS("def.longstartle")



#save(out_count, dt.window, observedata,bounds, pc_loading,pcr.ag, em.data,phenotype.dt,  file = "figuredataOctober.Rdata")
#########A  - karyotype modeling plot)
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
 ggplot(data = out_count[perm.status == "observed"], aes(x = group.general, y = N)) +
  geom_boxplot(width = 0.5, data = out_count[perm.status == "permutation"], 
               aes(x=group.general, y = N), outlier.colour = "NA") +
  geom_jitter(data = out_count[perm.status == "permutation"], 
              aes(x=group.general, y = N), alpha = 0.08, width = .2) +
  geom_point(shape = 23, 
             size = 5, 
             fill = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4","#F8766D", "#7CAE00", "#00B4C4", "magenta4","#F8766D", "#7CAE00", "#00B4C4", "magenta4","#F8766D", "#7CAE00", "#00B4C4", "magenta4") ) +
  facet_grid(cols = vars(inv), scales = "free")+
  theme_bw() + 
  
  xlab("Phenotypic Categories") +
  ylab( "# of Sig. Phenotypes (< 0.05) ")+
  theme(
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 10),
    #axis.title=element_text(size=12,face="bold"),
    axis.text.x = element_text( size = 12),
    axis.text.y = element_text(size = 12, face = "bold")) -> panelA


#B make co.enrichment figure
#now start making graph showing or vs proportion

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
  scale_x_continuous(breaks = seq(0,3, by = 0.5)) +
  theme(
    #legend.position = c(0.875, 2.05),
    legend.position = "none",
    legend.key.size = unit(.4, "cm") , 
    legend.text = element_text(size = 8),
    axis.text.x = element_text( size = 7, face = "bold"),
  ) -> panelB

####Part C- create sliding window figure
#create boundaries of our peaks
final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.2", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.2, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


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
             size = 3.0) +
  theme_bw()+
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos") -> panelC

### part D: make phenotype pca


lim <- 1.8
ggplot(data=pc_loading) +
  coord_equal() +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(data=pc_loading[Dim.1<0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=1) +
  geom_text(data=pc_loading[Dim.1>0], aes(x=x, y=y, label=pheno), size = 2.75, vjust=1, hjust=0) +
  
  geom_segment(aes(x=x, y=y, xend=Dim.1, yend=Dim.2), color="grey", alpha=0.5) +
  xlim(-lim,lim) + ylim(-lim,lim) +
  theme_bw() + xlab("PCA Dimension 1 Loadings (12.5%)") + ylab("PCA Dimension 2 Loadings (7.8%)") -> panelD

####  part E/F make pca loading figures
pcr.ag %>% 
  dplyr::select(V2,  mu1, se1, gt.name) %>%
  mutate(pc = 1)->
  ms1
names(ms1)[2:3] = c("mu", "s")

pcr.ag %>% 
  dplyr::select(V2,  mu2, se2, gt.name) %>%
  mutate(pc = 2) ->
  ms2
names(ms2)[2:3] = c("mu", "s")

rbind(ms1, ms2) %>%
 ggplot(data=., aes(x=as.factor(gt.name), y=mu)) +
  geom_errorbar(aes(x=as.factor(gt.name), ymin=mu-1.96*s, ymax=mu+1.96*s), width=.1) +
  geom_point() +
  theme_bw()  +
  xlab(NULL) +
  ylab("PCA Dimension 1\nProjection") +
  facet_wrap(~pc, ncol = 1) -> panelE

### part G - make phenotype validation figure

#part 1- emtrends startle decay

ggplot(em.data, aes(x = f1.background, y = scaled.minute.trend, ymin = asymp.LCL, ymax = asymp.UCL, color = inversion.st)) +
  geom_errorbar(width = 0.1, position=position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5)) +
  #facet_grid(.~f1.background, scales = "free_y")+
  theme_bw() + 
  #scale_color_manual(values = group.colors) +
  facet_grid("model.slope" ~ . ) +
  ylab("activity.slope") -> panelF
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
  theme_bw() -> panelG


panelA
panelB
panelC
panelD
panelE
panelF
panelG

joint.fig = (panelA+panelB+panelC)/(panelD+panelE+(panelF/panelG))

###make patchwork figures
ggsave(joint.fig, file =  "fig.d5.pdf", width = 14, height = 7)
#gtotal = (g1/ g2 / g3) | g4
#gtotal = gtotal +  plot_annotation(tag_levels = 'A')
#ggsave(gtotal,  file =  "fig.d6.png", width = 10, height = 7)


##plort deficiencies
data.frame(
label=c("k2a","k2b", "k5a","k5b","k9a"),
start=c(2175620,2242285,5147258,5073453,8958155),
end=c(2450829,2374023,5305636,5145500,9581740),
yaxis=1:5
) %>%
  ggplot(aes(
    xmin=start,
    xmax=end,
    ymin=yaxis-0.05,
    ymax=yaxis+0.05,
    fill=label
  )) +
  geom_rect() +
  xlim(0,20e6)
  