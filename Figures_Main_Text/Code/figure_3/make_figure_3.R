##### Prepare Figure 3
##### For this figure, i envision a 4 panel figure.

#### Panel 1: binom
#### Panel 2: wZa
#### Panel 3: LD
#### Panel 4: Phenotype
#### 

### load libraries
library(patchwork)
library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(data.table)
library(reshape2)
library(car)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(viridis)
library(foreach)
library(doMC)
library(lubridate)
library(forcats)
library(viridis)
registerDoMC(2)



##### Figure 1
### P.value distributions
output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)

#### load windows
#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
load(output_results_window)
###
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0 & chr.x == "2L",
                       list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
                       list(mod, chr.x=chr.x, locality, win.i, start, end)] %>%
  filter(locality == "VA_ch",
         mod=="aveTemp+year_factor")


win.out %<>% filter(locality == "VA_ch", mod=="aveTemp+year_factor", chr.x == "2L")
###win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

###win.minp.ag.ag <- win.minp.ag[perm!=0, 
## list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]

### basic MH plot
### Binomial p-value. 
mh.plot.binom <-
  ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="dashed") +
  #geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  #geom_line(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
  #          aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p)), size=0.2) +
  geom_ribbon(data=win.minp.ag,
              aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.45) +
  geom_point(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
             aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                 color=(rnp.pr)), size=1.1) +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab(expression(paste(-log[10], "(Window P)"))) 

ggsave(mh.plot.binom, file="mh.plot.binom.pdf", h=2.9, w=5)
####


mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="dashed") +
  #geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_line(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
            aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rnp.wZa.p)), size=0.8, color = "steelblue") +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab(expression(paste(-log[10], "(wZa P)"))) 

ggsave(mh.plot.wza, file="mh.plot.wza.pdf", h=2.9, w=5)

#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 
#### LD plot
load("../7.LD/merged.ld.Rdata")
### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B))

## bring in the inversion data -- for labeling SNPS.
inv_markers_id <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)

## Add inv/outlier label
ld_df %<>%
  mutate(A_status = case_when(SNP_A %in% inv_markers_id$V1 ~ "inv",
                              !(SNP_A %in% inv_markers_id$V1) ~ "out"),
         B_status = case_when(SNP_B %in% inv_markers_id$V1 ~ "inv",
                              !(SNP_B %in% inv_markers_id$V1) ~ "out"))  %>% 
  mutate(Comp_stat = paste(A_status, B_status, sep = "_")) %>%
  mutate(Comp_simplified = case_when(A_status == "inv" | B_status == "inv" ~ "vs_inv",
                                     A_status != "inv" & B_status != "inv" ~ "vs_out"))
############## triangle triangle triangle triangle triangle triangle triangle 
#### Plot LD triangle
ld_df %>%
  filter(R2 > 0.1,
         BP_A < BP_B) -> all_dat

ld_df %>%
  filter(R2 > 0.6,
         BP_A < BP_B) -> high_dat
############## triangle triangle triangle triangle triangle triangle triangle 
############## triangle triangle triangle triangle triangle triangle triangle 
############## triangle triangle triangle triangle triangle triangle triangle 
ggplot() +
  #geom_point(data=all_dat,
  #           aes(x=BP_A,
  #               y=BP_B,
  #               color = R2),
  #           size = 1, alpha = 0.5, shape = 15) +
  geom_point(data=high_dat,
             aes(x=BP_A,
                 y=BP_B,
                 color = R2),
             size = 0.09, alpha = 0.1, shape = 15) +
  theme_classic() +
  scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.80) ->
  ld_triag_plot

ggsave(ld_triag_plot, file = "raster_plot.png",
       width = 5,
       height = 4
       )
############## triangle triangle triangle triangle triangle triangle triangle 
############## triangle triangle triangle triangle triangle triangle triangle 
############## triangle triangle triangle triangle triangle triangle triangle 

#### Describe LD outliers

LD_GLM_outs <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

LD_GLM_outs %<>%
  mutate(tidy_annot = case_when(
    Annotation %in% c("3_prime_UTR_variant", "5_prime_UTR_variant", "non_coding_transcript_exon_variant") ~ "UTR",
    Annotation %in% c("downstream_gene_variant", "intergenic_region", "upstream_gene_variant") ~ "intergenic",
    Annotation %in% c("intron_variant") ~ "intron",
    Annotation %in% c("missense_variant") ~ "NS",
    Annotation %in% c("missense_variant&splice_region_variant", 
                      "splice_region_variant", 
                      "splice_region_variant&intron_variant") ~ "Splice",
    Annotation %in% c("synonymous_variant") ~ "S"))


LD_GLM_outs %>%
  group_by(tidy_annot ) %>%
  summarize(N= n()) %>%
  mutate(Nprop = prop.table(N)*100)

LD_GLM_outs %>%
  group_by(Transcript_biotype ) %>%
  summarize(N= n()) %>%
  mutate(Nprop = prop.table(N)*100)

LD_GLM_outs %>%
  filter(Transcript_biotype == "protein_coding") %>%
  group_by(Gene_Name ) %>%
  summarize(N= n()) %>%
  mutate(Nprop = prop.table(N)*100)


### pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno
### pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno
### pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno
### pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno
### pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno pheno

inv.dt <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/InversionsMap_hglft_v6_inv_startStop.txt")

setnames(inv.dt, "chrom", "chr")


load("/project/berglandlab/alan/all_gwas_glm_win3.Rdata")

#gwas.win.o[,grm:=!grepl("nogrms", gwas.pheno)]
#gwas.win.o.ag <- gwas.win.o[,list(n=sum(pa<.05 & or>1)), list(chr, start, end, grm, glm.perm.n, inv, win.i)]
gwas.win.o.ag.ag <- gwas.win.o.ag[,list(n=mean(n), sd=sd(n)), list(chr, start, end, grm, perm=I(0!=glm.perm.n), inv)]


ggplot() +
  #geom_rect(data=inv.dt[which(inv.dt$invName == "2Lt"),],
  #          aes(xmin=start/1e6, xmax=stop/1e6, ymin=-1, ymax=20), color="grey", alpha=.5) +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_ribbon(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
              aes(x=I(start), 
                  ymax=n+2*sd,
                  ymin=0), 
              color="black", alpha = 0.6) +
  geom_point(data=gwas.win.o.ag.ag[perm==F][grm==F][chr=="2L"],
             #aes(x=I(start/2+end/2)/1e6, y=n), size = 1, color ="red") +
             aes(x=I(start), y=n), size = 1.7, color ="red", alpha = 0.5) +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
            #aes(x=I(start/2+end/2)/1e6, y=n), color="black") +
            aes(x=I(start), y=n), color="black") +
  xlim(0,23500000) +
  facet_grid(grm~chr) +
  theme_bw() -> pheno_plot

ggsave(pheno_plot, file = "pheno_plot.pdf", h=2.9, w=5)

unique(gwas.win.o[glm.perm.n==0][chr=="2L"][pa<.0005][inv==T]$gwas.pheno)
unique(gwas.win.o[glm.perm.n==0][chr=="2L"][inv==T]$gwas.pheno)


#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 
#### Tajima's D plus wZA
#### Tajima's D plus wZA
#### Tajima's D plus wZA
#### Tajima's D plus wZA
#### Tajima's D plus wZA
#### Tajima's D plus wZA
#### Tajima's D plus wZA

wza_dat <- win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0]
wza_dat %<>%
  mutate(wza=-1*log10(rnp.wZa.p)) %>%
  dplyr::select(BIN_START=start, wza) %>%
  melt(id = "BIN_START") %>%
  mutate(pop = "cm",
         type = "pool")

load("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/pi_D_datforplot.Rdata")
D_dat = sub_pi_d_parsed %>%
  filter(resolution == "W_100000",
         variable == "TajimaD") %>% 
  dplyr::select(BIN_START, pop, type, variable, value)

rbind(wza_dat, D_dat) -> wza_and_D

wza_and_D %>%
  ggplot(
    aes(
      x=BIN_START,
      y=value,
      color=pop,
      linetype=type
    )) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_line() +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  wza_plus_D_plot

ggsave(wza_plus_D_plot, file = "wza_plus_D_plot.pdf")


#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 
#  The multipopulation test
#### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### 

#load("/project/berglandlab/alan/joint_window.Rdata")
#
##inv.dt <- fread("Overwintering_18_19/Inversions/InversionsMap_hglft_v6_inv_startStop.txt")
#inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions#/InversionsMap_hglft_v6_inv_startStop.txt"
#### load suppl data
#inv.dt <- fread(inversion_map)
#
#
#setnames(inv.dt, "chrom", "chr")
#
#
#win.out = win.out[chr.x == "2L"]
#
#win.out.ag <- win.out[,list(n=sum(pa<.1), nl=sum(fet.p<.0005), pops=paste(unique(locality.y[fet.p#<.0005]), collapse="_")), list(thr, win.i, perm, start.bp, stop.bp, chr=chr.x)]
#win.out.ag.ag <- win.out.ag[,list(n=mean(nl), lci=quantile(nl, .025), uci=quantile(nl, .975), pops#=pops[1]), list(thr, win.i, perm=I(perm!=0), start.bp, stop.bp, chr)]
#
#g1 <-
#  ggplot() +
#  geom_vline(data=inv.dt[chr == "2L"], aes(xintercept=start, linetype=invName)) +
#  geom_vline(data=inv.dt[chr == "2L"], aes(xintercept=stop, linetype=invName)) +
#  geom_ribbon(data=win.out.ag.ag[perm==T][thr==0.05], aes(x=start.bp/2 + stop.bp/2, ymin=lci, ymax#=uci), fill="black", alpha=.5) +
#  geom_point(data=win.out.ag.ag[perm==F][thr==0.05], aes(x=start.bp/2 + stop.bp/2, y=n, color=pops#), size = 2, alpha = 0.5) +
#  facet_grid(~chr) +
#  theme_bw()
#
#ggsave(g1, file="./multipop_overlap.pdf", h=2.9, w=5)
#

