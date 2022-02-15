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

#### Find names of SNPs for the lower triangle of LD




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
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start/1e6), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop/1e6), linetype="dashed") +
  geom_ribbon(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
              aes(x=I(start/2+end/2)/1e6, 
                  ymax=n+2*sd,
                  ymin=0), 
              color="black", alpha = 0.6) +
  geom_point(data=gwas.win.o.ag.ag[perm==F][grm==F][chr=="2L"],
             aes(x=I(start/2+end/2)/1e6, y=n), size = 1, color ="red") +
  geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
            aes(x=I(start/2+end/2)/1e6, y=n), color="black") +
  #geom_line(data=gwas.win.o.ag.ag[perm==T][grm==F][chr=="2L"],
  #      aes(x=I(start/2+end/2)/1e6, y=n-2*sd), color="black", linetype="dashed") +
  facet_grid(grm~chr) +
  theme_bw() -> pheno_plot

ggsave(pheno_plot, file = "pheno_plot.pdf", h=2.9, w=5)

unique(gwas.win.o[glm.perm.n==0][chr=="2L"][pa<.0005][inv==T]$gwas.pheno)
unique(gwas.win.o[glm.perm.n==0][chr=="2L"][inv==T]$gwas.pheno)



