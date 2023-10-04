### New Figure Sliding window Phenos!
### 
library(tidyverse)
library(data.table)
library(FactoMineR)
library(ggrepel)
library(reshape2)
library(gmodels)

load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/FIGURES_TABLES_FILES/pub.Figures_Supp_Text/Revision_Figures/gwas_ld_prune_object.Rdata")

invs <- fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/FIGURES_TABLES_FILES/pub.Figures_Supp_Text/Revision_Figures/InversionsMap_hglft_v6_inv_startStop.txt")
names(invs)[1] = "chr"

  ggplot() + 
  geom_rect(
  data=invs,
  aes(
    xmin=start/1e6,
    xmax=stop/1e6,
    ymin = 0,
    ymax = 14,
    fill=invName
  ), alpha = 0.2) +
  geom_ribbon(
    data=dt.window,
    aes(
    x=(start+end)/2e6,
    ymax=upperbound,
    ymin = 0), alpha = 0.4) + 
  geom_point(
    data=dt.window,
    aes(
    x=(start+end)/2e6,
    y=avg.N,
    color=`sig.v.per`
  ))  +
    theme_bw() +
    theme(legend.position = "top") +
  facet_grid(~chr, scales = "free_x")

###
# newobserved = observedata
  # newobserved = observedata

  load("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/FIGURES_TABLES_FILES/pub.Figures_Supp_Text/Revision_Figures/ldprunecrescent")
  
  bounds$inv = gsub(FALSE, "Not-Inverted", bounds$inv)
  bounds$inv = gsub(TRUE, "Inverted", bounds$inv)
  observedata$inv = gsub(FALSE, "Not-Inverted", observedata$inv)
  observedata$inv = gsub(TRUE, "Inverted", observedata$inv)
  #bounds = bounds[chr == "2L"]
  #observedata = observedata[chr == "2L"]
  
  g1 =ggplot() +
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
    # scale_x_continuous(breaks = seq(0,3, by = 0.5)) +
    xlim(0,2.5) +
    theme(
      #legend.position = c(0.875, 2.05),
      legend.position = "none",
      axis.title.x=element_blank(),
      legend.key.size = unit(.4, "cm") ,
      legend.text = element_text(size = 8),
      axis.text.x = element_text( size = 7, face = "bold"),
    ) +
    ggtitle("NoGRM")
  
###
### --- aditions to figure 7!
  bounds$inv = gsub(FALSE, "Not-Inverted", bounds$inv)
  bounds$inv = gsub(TRUE, "Inverted", bounds$inv)
  observedata$inv = gsub(FALSE, "Not-Inverted", observedata$inv)
  observedata$inv = gsub(TRUE, "Inverted", observedata$inv)
  bounds = bounds[chr == "2L"]
  observedata = observedata[chr == "2L"]

  g2 =ggplot() +
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
    # scale_x_continuous(breaks = seq(0,3, by = 0.5)) +
    xlim(0,2.5) +
    theme(
      #legend.position = c(0.875, 2.05),
      legend.position = "none",
      axis.title.x=element_blank(),
      legend.key.size = unit(.4, "cm") ,
      legend.text = element_text(size = 8),
      axis.text.x = element_text( size = 7, face = "bold"),
    ) +
    ggtitle("NoGRM")
##### .
##### 
##### 
##### 
  final.windows.pos = 
    data.frame(win.name = c(#"left", 
      "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6"
      #, "right" 
    ),
    mid = c(#2.2, 
      3.1, 4.7, 5.2, 6.1, 6.8 , 9.6
      #, 13.1
    ),
    chr = "2L"
    ) %>%
    mutate(start = (mid-0.2)*1e6 ,
           end  = (mid+0.2)*1e6  )
  
  
  ggplot() +
    geom_rect(data = final.windows.pos, aes(xmin=start/1e6 , xmax =end/1e6, ymin = 0, ymax = 14),
              fill = "gold", alpha = 0.6) +
    geom_vline(xintercept = 2225744/1e6) +
    geom_vline(xintercept = 13154180/1e6) +
      geom_ribbon(
      data=dt.window[chr=="2L"],
      aes(
        x=(start+end)/2e6,
        ymax=upperbound,
        ymin = 0), alpha = 0.4) + 
    geom_point(
      data=dt.window[chr=="2L"],
      size = 2.5,
      aes(
        x=(start+end)/2e6,
        y=avg.N,
        color=`sig.v.per`
      ))  +
    scale_color_manual(values = c("purple","darkcyan")) +
    theme_bw() +
    theme(legend.position = "top") +
    facet_grid(~chr, scales = "free_x")
  
  
  ### Panel -PCA
  dt.window %>% 
    filter(chr == "2L") %>%
    filter(sig.v.per == TRUE) %>%
    filter(start > 2225744 & end < 13154180) %>%
    .$gwas.pheno -> Phenotab.sw

###
###
sig.phenos = fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/FIGURES_TABLES_FILES/pub.Figures_Main_Text/revision_genetics/sigPhenos.txt")
phenos = readRDS("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/CODE/13.Phenotype_Analysis/DATA/wideform.fixed.phenotable.RDS")
as.data.frame(phenos) -> phenos
cols = which(names(phenos) %in% sig.phenos$phenoname)

phenos[,c(cols)] -> dat.for.PCA
row.names(dat.for.PCA) = phenos$ral_id

dat.for.PCA %>%
  PCA(grap = F) ->
  PCA.updated

PCA.updated$var$coord %>%
  as.data.frame() %>% 
  mutate(pheno = row.names(.)) %>%
  separate(pheno, into = c("pheno","etc"), sep = "_") -> dat.for.plot
setDT(dat.for.plot)

  ggplot()+
  geom_segment(
    data=dat.for.plot,
    aes(x=0, y=0, xend=Dim.1, yend=Dim.2),
    arrow=arrow(length=unit(0.05,"cm"))
    ,color="grey", alpha = 1.0, linewidth = 0.2) +
    geom_text(data=dat.for.plot[Dim.1<0], aes(x=Dim.1, y=Dim.2, label=pheno), size = 2.5, vjust=1, hjust=1) +
    geom_text(data=dat.for.plot[Dim.1>0], aes(x=Dim.1, y=Dim.2, label=pheno), size = 2.5, vjust=1, hjust=0) +
    xlim(-0.55,0.55) +
  theme_bw()

## 
  dat.for.plot %>%
    filter( Dim.1 < 0 & Dim.2 < 0 ) %>%
    .$pheno %>% table %>% as.data.frame()

  dat.for.plot %>%
    filter( Dim.1 < 0 & Dim.2 > 0 ) %>%
    .$pheno %>% table %>% as.data.frame()
  
  dat.for.plot %>%
    filter( Dim.1 > 0 & Dim.2 > 0 ) %>%
    .$pheno %>% table %>% as.data.frame()

  dat.for.plot %>%
    filter( Dim.1 > 0 & Dim.2 < 0 ) %>%
    .$pheno %>% table %>% as.data.frame()

######. Corrlations
inv.sts <- fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/FIGURES_TABLES_FILES/pub.Figures_Main_Text/revision_genetics/dgrp.inverion.status.txt")
######
  PCA.updated$ind$coord %>%
    as.data.frame() %>%
    mutate(DGRP_line = 
             paste("DGRP", row.names(.), sep = "_")) %>%
  left_join(inv.sts) %>%
  filter(inv2lt %in% c("ST","INV")) ->
  mapped.data.toinvs
  setDT(mapped.data.toinvs)
mapped.data.toinvs %>%
  dplyr::select(inv2lt, DGRP_line, Dim.1, Dim.2 ) %>%
  reshape2::melt(id = c("inv2lt", "DGRP_line")) %>%
  group_by(inv2lt, variable) %>%
  summarise(mean = ci(as.numeric(value))[1],
            low = ci(as.numeric(value))[2],
            high = ci(as.numeric(value))[3],
            ) %>%
  ggplot(aes(
    x=inv2lt,
    y=mean,
    ymin=low,
    ymax=high
  )) + geom_errorbar(width = 0.2) +
  geom_point() + 
  facet_grid(variable~.) +
  theme_bw()

  
  t.test(mapped.data.toinvs[inv2lt == "ST"]$Dim.1,
         mapped.data.toinvs[inv2lt == "INV"]$Dim.1)

  t.test(mapped.data.toinvs[inv2lt == "ST"]$Dim.2,
         mapped.data.toinvs[inv2lt == "INV"]$Dim.2)
  
  t.test(mapped.data.toinvs[inv2lt == "ST"]$Dim.3,
         mapped.data.toinvs[inv2lt == "INV"]$Dim.3)
  
  