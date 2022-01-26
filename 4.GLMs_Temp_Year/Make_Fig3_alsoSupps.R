### Prepare Panels for figure 4. 
### 

library(tidyverse)
library(vroom)
library(tidyverse)
library(data.table)
library(tidyr)
library(viridis)
library(patchwork)
library(magrittr)
library(RColorBrewer)
library(rstatix)


#### Part 1 --- Describe the distribution of p.values

### ---> SEE special script ---> 

##### Temperature Figure -- Not part of figure 3

#github_addr <- "/Users/jcbnunez/Documents/GitHub"
### Part 1: load data
#load(paste(github_addr,"/Cville-Seasonality-2016-2019/4.GLMs_Temp_Year/temperatureAverage_yearFactor_GLM/weatherAve#.Rdata", sep =""))
#
#metadat <- vroom(paste(github_addr,"/Cville-Seasonality-2016-2019/4.GLMs_Temp_Year/temperatureAverage_yearFactor_GLM#/DEST_10Mar2021_POP_metadata.csv", sep =""))
#
#names(weather.ave)[1] = "sampleId"
#
#weather.ave %>%
#  left_join(metadat) %>%
#  mutate(date_posix = as.Date(collectionDate, format = "%m/%d/%Y")) %>%
#  ggplot(aes(x = date_posix, 
#             y = aveTemp/10,
#             group = year,
#             fill = aveTemp/10
#             )) +
#  geom_line() +
#  geom_point(shape = 21, size = 1.5) +
#  scale_x_date(date_labels = "%y") +
#  scale_fill_gradient2(low = "steelblue", high = "firebrick", midpoint = 15) +
#  theme_bw() +
#  theme(legend.position = "top") +
#  facet_grid(.~locality, scales = "free_x",  space = "free_x")
#
#########
######### Plot the giant manhattan plot
#

output_results_window <- "/project/berglandlab/thermal_glm_dest/window_analysis_output.nested.qb.Rdata"
inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

### load suppl data
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

ace.dt <- data.table(chr.x="3R", pos=13174333/2  + 13274333/2)


#### load windows
#load("/Users/alanbergland/Documents/GitHub/Overwintering_18_19/temperatureAverage_yearFactor_GLM/3.windowAnalysis/window_analysis_output.Rdata")
# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.Rdata ~/.
#load("~/window_analysis_output.nested.Rdata")

#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/window_analysis_output.nested.qb.Rdata ~/.
load(output_results_window)

### get sig thresholds from perms
## win.minp.ag <- win.out[pr==0.01 & nSNPs>100,
##         list(q05=quantile(wZa.p, 0.025, na.rm=T), q5=quantile(rbinom.p, .05, na.rm=T), .N),
##         list(perm, mod, locality)]$N %>% mean
## win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

## win.minp.ag.ag <- win.minp.ag[perm!=0, list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]

### v2 -- stablish significance
###
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0,
                       list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
                       list(mod, chr.x=chr.x, locality, win.i, start, end)]

###win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

###win.minp.ag.ag <- win.minp.ag[perm!=0, 
## list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]


### Part 1.  ---> the manhattan plot

### basic MH plot
mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_ribbon(data=win.minp.ag[mod=="aveTemp+year_factor"],
              aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.75) +
  geom_point(data=win.out[mod=="aveTemp+year_factor"][pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
             aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                 color=(rnp.pr)), size=.95) +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(locality~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.box = "vertical",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab("-log10(Window p)")

ggsave(mh.plot.wza, file="nested_qb.png", h=8.5, w=11)

