#setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/REPRODUCE_Fig6")
###Inversion paper figure
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(corrplot)
library(patchwork)
library(ggrepel)
#install.packages('doParallel', repos='http://cran.us.r-project.org')
#use doParallel package to register multiple cores that can be used to run loops in parallel
registerDoParallel(10)
#setwd("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/REPRODUCE_Fig6")
#load in data 
#broad coenrichment data
load("figuredatajuly.Rdata")
#karyotype modeling data
#out_count = readRDS("karyotype.modeling.RDS")
#coloc heatmap data
#alldata = readRDS( "coloc.permutations.RDS")
#setwd("/Users/supad/OneDrive/Documents/Bergland Research/R_data_objects/July_2022_objects/")
#slding window data
# groupby.perm = readRDS("newmodelplotdata.RDS")
#dt.window = readRDS("slid.window.2.5.plot.RDS")

#phenotype figure data
# groupbysnp = readRDS( "hot.cold.phenoeffect.RDS")
#groupbysnp = readRDS("hot.cold.phenoeffect.RDS")

#save(out_count, alldata, groupbysnp, dt.window, observedata,bounds, file = "figuredatajuly.Rdata")
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
#merge in beat count
#add in labels of how many permutations each observed beats.
beats = foreach(f = unique(out_count$chr), .combine = "rbind") %do% {
  # f = "2L"
  tmp = out_count[chr == f]
  group.o = foreach(g = unique(out_count$group.general), .combine ="rbind") %do% {
    # g = "B."
    tmp2 = tmp[group.general == g]
    #show threshold permutations have to beat
    threshold = tmp2[perm.status == "observed"]$N
    beatcount = (tmp2[perm.status == "permutation"][N > threshold])
    beatcount = dim(beatcount)[1]
    o = data.table(
      chr = f,
      group.general = g,
      beatcount = beatcount
    )
    o
  }
  group.o
}
merg.inv.data = merge(out_count, beats, by = c("chr", "group.general"))
#now show number of perms beaten by oversved = 100 - perms beating it
merg.inv.data$beatcount = 100 - merg.inv.data$beatcount
head(out_count)

#try graph again
g1 =  ggplot(data = merg.inv.data[perm.status == "observed"], aes(x = group.general, y = N, label = beatcount)) +
  geom_boxplot(width = 0.5, data = merg.inv.data[perm.status == "permutation"], 
               aes(x=group.general, y = N), outlier.colour = "NA") +
  geom_jitter(data = merg.inv.data[perm.status == "permutation"], 
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
    axis.title=element_text(size=12,face="bold"),
    axis.text.x = element_text( size = 12),
    axis.text.y = element_text(size = 12, face = "bold")) +
  geom_text_repel(nudge_y = 2)


#B make co.enrichment figure
#now start making graph showing or vs proportion

bounds$inv = gsub(FALSE, "Not-Inverted", bounds$inv)
bounds$inv = gsub(TRUE, "Inverted", bounds$inv)
observedata$inv = gsub(FALSE, "Not-Inverted", observedata$inv)
observedata$inv = gsub(TRUE, "Inverted", observedata$inv)

g2 = ggplot() +
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
  ) 
#B make co.enrichment figure


####Part C- create sliding window figure
#create boundaries of our peaks
final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


g3 = ggplot() +
  #geom_rect(data=outlier_haplowins, aes(xmin=start/1e6, xmax=end/1e6, ymin=-1, ymax=100),fill ="lightgoldenrod1", alpha=.5) +
  geom_rect(data = final.windows.pos,
            aes(xmin=start/1e6, xmax = end/1e6,
                ymin = 0, ymax = 85), 
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
             size = 1.2) +
  scale_color_manual(values = c("slateblue4", "violet")) +
  theme_bw()+
  # ylim(0,25) +
  xlab("Position on Chromosome") +
  ylab( "# of Co-enriched Phenos")

## how many phenotypes in Msp300
dt.window[perm.st == "observed"] %>%
  arrange(-avg.N)

####D- make phenotype effect figure
#reorder peak groups for the facet
groupbysnp$peak = factor(groupbysnp$peak, 
                         levels = c("inv.start", "win_4.6", "win_5.1", "win_6.8", "win_9.5", "inv.stop"))

section = groupbysnp[peak == "win_5.1"][phenotype %in% c(
  "ActivityLevel_Standard-BasalActivity_F",
  "NegativeGeotaxis_MSB_Treatment_female",
  "StartleResponse_standard_female",
  "StarvationResistance_standard_male",
  "MeanElutionTime_DevelopmentofTolerance_female",
  "FreeGlycerolLevels_LowGlucoseDiet_male"
)]

ggsave(g3, file = "g3.pdf", h = 3.5, w =12)

g4 = ggplot(data = section, aes(x = as.factor(temp), y = mean.deviation, group = phenotype, color = group.general, label = label)) +
  geom_line(size = 1)+
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00B4C4", "magenta4"))  +
  xlab("temperature when this variant is at higher frequency") +
  ylab("Distance from the mean") +
  facet_wrap(~phenotype)+
  theme(
    legend.position = "none"
  ) 

####plot D part 2= a table of 
#### Emake coloc heatmap figure
peaks = unique(alldata$peak)
outpeaks = foreach(f = peaks) %do% {
  f = peaks[4]
  peak.interest = f
  print(peak.interest)
  #filter to data from a certain peak
  peakdata = alldata[peak == peak.interest]
  peakdata = peakdata[, -1]
  #use dcast to create wideform-h4 data
  widedata = subset(peakdata, select=c(meanh4, pheno1.fixed, pheno2.fixed))
  widedata = dcast(widedata, pheno1.fixed ~ pheno2.fixed, value.var = "meanh4")
  rows = widedata[[1]]
  widedata = widedata[,-1]
  M = as.matrix(sapply(widedata, as.numeric))
  
  rownames(M) = rows
  #do the same to create a matrix of p values
  pdata =  subset(peakdata, select=c(pvalue, pheno1.fixed, pheno2.fixed))
  pdata = dcast(pdata, pheno1.fixed ~ pheno2.fixed, value.var = "pvalue")
  rows = pdata[[1]]
  pdata = pdata[,-1]
  P = as.matrix(sapply(pdata, as.numeric))
  rownames(P) = rows
  #plot a heatmap using colors, ordered by aod
  filename = paste0("sig.contrasts", peak.interest, ".png")
  png(filename, width = 480, height = 480)
  #  heatmap(M,  symm = T, cexRow = 0.3, cexCol = 0.3)
  
  corrplot(M, type = "lower" , method = "color",order = "AOE",  
           tl.cex = 0.55, tl.col = "black", tl.offset = 1,
           p.mat = P, sig.level = 0.5, insig = "pch", pch.col = "black", pch.cex = 0.5)
  dev.off()
}
x = c(rep("red", 20), rep("blue", 19))
###make patchwork figures
ggsave((p1 / p2) | (p3 / p4), file =  "fig.d4.png", width = 14, height = 7)
gtotal = (g1/ g2 / g3) | g4
gtotal = gtotal +  plot_annotation(tag_levels = 'A')
ggsave(gtotal,  file =  "fig.d6.png", width = 10, height = 7)
