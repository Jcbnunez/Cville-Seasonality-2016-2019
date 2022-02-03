### Make Figure 3
### 
### For figure 3 I envision a 4xn panel figure. The the first row if panels will showcase the temperature data
### The point of this panel is to show the temperature arcs
### The second panel will show the birds eye view of the distribution of LRT p-values for the 7 pops
### The third line will show the enrrichment analysis for charlottesville
### the forth line will show the manhattan plot for charlottesville 
### 

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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### This is built using "getEnvironmentMatrix.aveTemp_yearFactor.R"
### SUPPLEMENT NOW
load("./weatherAve.Rdata")
names(weather.ave)[1] = "sampleId"
samps <- fread("/project/berglandlab/DEST_Charlottesville_TYS/DEST_metadata/DEST_10Mar2021_POP_metadata.csv")
samps <- samps[set!="dgn"]
samps$locality[grep("UA_od", samps$locality )] = "UA_Ode"

left_join(weather.ave[,-3], samps, by = "sampleId" ) -> weather_samps

weather_samps %>%
  mutate(date_posix = as.Date(collectionDate, format = "%m/%d/%Y")) %>%
  ggplot(aes(x = date_posix, 
             y = aveTemp/10,
             group = year,
             fill = aveTemp/10
  )) +
  geom_line() +
  geom_point(size = 1.2, shape = 21, color = "black") +
  scale_fill_gradient2(low = "steelblue", high = "firebrick", midpoint = 15) +
  theme_bw() +
  scale_x_date(date_labels = "%y") +
  theme(legend.position = "top") +
  facet_grid(.~locality, scales = "free_x") ->
  city_temp

ggsave(city_temp, 
       file = "city_temp.pdf",
       width = 9,
       height = 2.2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### General GLM results
### SUPPLEMENT NOW

files_to_read <- system("ls ./p_val_extractor_w_CRM | grep '.Rdata' ", intern = T)


collect_counts = list()
collect_enrrich = list()
for(i in 1:length(files_to_read)){
  print(files_to_read[i])
  
  load( paste("/scratch/yey2sn/Overwintering_ms/4.GML_plots/p_val_extractor_w_CRM/",
              files_to_read[i],
              sep = ""))
  
  collect_counts[[i]] = count_df
  collect_enrrich[[i]] = enrrich_df
  
}

#quick post-processing
counts_df = do.call(rbind, collect_counts)
counts_df %<>%
  mutate(proportion = Tresh/N)

counts_df %<>%
  mutate(fill_for_plot = paste(inversion_pos,CRM_pos, sep = "")) 
counts_df$fill_for_plot = gsub("NA", "", counts_df$fill_for_plot )
counts_df$fill_for_plot = gsub("^$", "All", counts_df$fill_for_plot )
counts_df$fill_for_plot %>% table

counts_df$chr[grep("genome", counts_df$chr)] = "G"

dat_in = counts_df %>%
  filter(#pop == "VA_ch",
    p_tresh == 0.05,
    fill_for_plot %in% c("All", "inv", "non.inv")) 

dat_in$fill_for_plot[which(dat_in$fill_for_plot == "All")] = "3.All"
dat_in$fill_for_plot[which(dat_in$fill_for_plot == "inv")] = "1.inv"
dat_in$fill_for_plot[which(dat_in$fill_for_plot == "non.inv")] = "2.non.inv"

ggplot() +
  geom_violin(data=dat_in[which(dat_in$category == "Per"),],
              aes(x=chr,
                  y=proportion,
                  fill = fill_for_plot),
              alpha = 0.7
  ) +
  geom_point(data=dat_in[which(dat_in$category == "Obs"),],
             aes(x=chr,
                 y=proportion,
                 fill = fill_for_plot),
             size = 2.3,
             shape = 23,
             color = "black",
             position = position_dodge2(w = 0.95)) +
  facet_grid(~pop, scales = "free_x", space = "free") +
  ylab(expression( paste("%", italic(P[LRT]<0.05)) ) ) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "top")-> 
  p_tresh_plot

ggsave(p_tresh_plot, 
       file = "p_tresh_plot.pdf",
       width = 9,
       height = 2.5)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### Manhattan plot
####

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
win.minp.ag <- win.out[pr==0.05 & nSNPs>100 & perm!=0,
                       list(lci=quantile(rbinom.p, 0.025, na.rm=T), uci=quantile(rbinom.p, .975, na.rm=T), .N),
                       list(mod, chr.x=chr.x, locality, win.i, start, end)] %>%
  filter(locality == "VA_ch",
         mod=="aveTemp+year_factor")

win.out %<>% filter(locality == "VA_ch", mod=="aveTemp+year_factor")
###win.minp.ag[locality=="VA_ch"][chr.x=="2L"][mod=="aveTemp+year_factor"]

###win.minp.ag.ag <- win.minp.ag[perm!=0, 
## list(q5.minp=quantile(minp, .05, na.rm=T), min.q5=min(q5, na.rm=T)), list(locality, mod, chr.x)]

### basic MH plot
### Binomial p-value. 
mh.plot.wza <-
  ggplot() +
  geom_vline(data=inv.dt, aes(xintercept=start/1e6, linetype=invName)) +
  geom_vline(data=inv.dt, aes(xintercept=stop/1e6, linetype=invName)) +
  #geom_vline(data=ace.dt, aes(xintercept=pos/1e6), color="orange", size=2, alpha=.75) +
  geom_ribbon(data=win.minp.ag,
              aes(x=(start/2 +end/2)/1e6, ymin=-1*log10(uci), ymax=-1*log10(lci)),
              color="grey", fill="grey", alpha=0.75) +
  geom_point(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0],
             aes(x=(start/2 +end/2)/1e6 , y=-1*log10(rbinom.p),
                 color=(rnp.pr)), size=.95) +
  geom_hline(yintercept=-log10(.01/1800)) +
  #geom_hline(data=win.minp.ag.ag[mod=="aveTemp+year_factor"], aes(yintercept=-log10(min.q5))) +
  facet_grid(~chr.x, scales="free") +
  theme_bw() +
  scale_color_viridis(option="G", direction = -1) +
  theme(axis.text.x=element_text(angle=0),
        legend.box = "vertical",
        legend.key.size=unit(1/8, 'in'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  labs(color="Prop. top 1%", linetype="Inversion") +
  xlab("Position (Mb)") +
  ylab(expression(paste(-log[10], "(Window P)"))) 

ggsave(mh.plot.wza, file="nested_qb.pdf", h=1.7, w=9)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Panel 
### VA GLM results enrriched

files_to_read <- system("ls ./P_value_VA_windows | grep '.Rdata' ", intern = T)


collect_counts = list()
collect_enrrich = list()
for(i in 1:length(files_to_read)){
  print(files_to_read[i])
  
  load( paste("/scratch/yey2sn/Overwintering_ms/4.GML_plots/P_value_VA_windows/",
              files_to_read[i],
              sep = ""))
  
  collect_counts[[i]] = count_df
  collect_enrrich[[i]] = enrrich_df
  
}

enrrich_df = do.call(rbind, collect_enrrich)


dat_in = enrrich_df %>% filter(!is.infinite(OR), 
                               !is.na(OR), 
                               p_tresh == 0.05, 
                               chr == "2L") %>%
  mutate(win.med = (start.bp+stop.bp)/2 ) %>% 
group_by(chr, category, Tidy_Consequence, win.med) %>%
  summarize(
            N = n(),
            OR_m = mean(OR),
            OR_sd = sd(OR),
            OR_l = quantile(OR, 0.05),
            OR_h = quantile(OR, 0.95))


ggplot() +
  geom_vline(xintercept = 2225744) + 
  geom_vline(xintercept = 13154180) + 
  geom_ribbon(data=dat_in[which(dat_in$category == "Per"),],
              aes(x=win.med,
                  ymax=log2(OR_h),
                  ymin=log2(OR_l),
                  #fill = inversion_pos
                  ), alpha = 0.7) +
  geom_point(data=dat_in[which(dat_in$category == "Obs"),],
             aes(x=win.med,
                 y=log2(OR_m),
                 #fill = inversion_pos
                 ),
             size = 1.1,
             shape = 23,
             color = "red",
             #position = position_dodge2(w = 0.75)
             ) +
  #facet_grid(~chr, scales = "free_x", space = "free") +
  ylab("OR") +
  #scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  geom_hline(yintercept = log2(1)) +
  theme(legend.position = "top") +
  facet_grid(Tidy_Consequence~., scales = "free_x") ->
  p_tresh_plot_enrrich

ggsave(p_tresh_plot_enrrich, 
       file = "p_tresh_plot_enrrich.png",
       width = 9,
       height = 9)

########## Functional Enrrichment
########## Functional Enrrichment
########## Functional Enrrichment
########## Functional Enrrichment
########## Functional Enrrichment
########## Functional Enrrichment
########## Functional Enrrichment

load("./VA_ch.0.05.2L.gene.enrrich_count.Rdata")

genes <- unique(Local_enrrichment_df$Symbol)

##dat_in <- Local_enrrichment_df 

annotation_file <- "/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata"
load(annotation_file)

### Make gene addresses
annotation_dmel_sp %>%
  group_by(chr, Symbol) %>%
  summarise(Median_pos = median(pos)) ->
  position_genes

### Generate Id list
prot_cod_ids <- unique(Local_enrrichment_df$Symbol)

Local_enrrichment_df %>%
  group_by(Symbol, category, ith, Analysis) %>%
  summarise(All_SNPs_tot = sum(All_SNPs),
            Outliers_tot = sum(Outliers)) %>% 
  mutate(proportion = Outliers_tot/All_SNPs_tot) ->
  Local_enrrichment_df_summarize

###
Local_enrrichment_df_summarize %>%   
  left_join(position_genes) -> 
  Local_enrrichment_df_summarize_pos

Local_enrrichment_df_summarize_pos %<>% 
  group_by(category, Symbol, Median_pos, Analysis, All_SNPs_tot) %>%
  summarise(Median = median(proportion),
            prop_l = quantile(proportion, 0.05),
            prop_h = quantile(proportion, 0.95),
            prop_sd = sd(proportion)
  ) 

Local_enrrichment_df_summarize_pos %>%
  filter(category == "Per") %>% 
  dplyr::select(Median_pos, category, Analysis, Symbol, top_9X=prop_h) ->
  top_9X

left_join(Local_enrrichment_df_summarize_pos, top_9X[,-2]) -> 
  Local_enrrichment_df_summarize_pos_top

dat_in = Local_enrrichment_df_summarize_pos_top %>%
  mutate(Cat_pval = 1-pbinom(Median*All_SNPs_tot, All_SNPs_tot, .05)) %>%
  filter(All_SNPs_tot >= 5) #%>%
#  mutate(Cat_pval_fdr = p.adjust(Cat_pval))
#.[complete.cases(.$Cat_pval),]

#save(dat_in, file = "./checkpoint_genes_of_interest.Rdata")
load("./checkpoint_genes_of_interest.Rdata")

dat_in$Analysis %>% table

dat_in %>% filter(Cat_pval < 1e-10,
                  Median>top_9X,
                  Analysis == "All") %>%
  .$Symbol %>% unique

dat_in %>% filter(Cat_pval < 1e-10,
                  Median>top_9X,
                  Analysis == "All",
                  Median_pos > 2225744,
                  Median_pos < 13154180) %>%
  .$Symbol %>% unique

dat_in %>% filter(Cat_pval < 1e-10,
                  Median>top_9X,
                  Analysis == "All",
                  Median_pos > 2225744,
                  Median_pos < 13154180) %>%
  .$Symbol %>% unique %>%
  .[grep("CG", ., invert = T)] %>%
  .[grep("lncRNA", ., invert = T)] 
  


library(ggrepel)
dat_in %>% filter(Analysis == "All")  -> dat_in_plot
ggplot() +
  #geom_point(
  #  data=dat_in[which(dat_in$category == "Per"),],
  #  aes(
  #  x=Median_pos,
  #  y=-1*log10(Cat_pval),
  #  #ymin = prop_l,
  #  #ymax = prop_h,
  #  fill = category,
  #  #color = category
  #), alpha = 0.3) +
  geom_point(
    data=dat_in_plot[which(dat_in_plot$category == "Obs"),],
    aes(
      x=Median_pos,
      y=-1*log10(Cat_pval),
      fill = category), 
    size = 0.7,
    shape = 23,
    alpha = 0.6
  )  +
  geom_point(
    data=dat_in_plot[which(dat_in_plot$category == "Obs" &
                             dat_in_plot$Median > dat_in_plot$top_9X &
                          dat_in_plot$Cat_pval < 1e-10),],
    aes(
      x=Median_pos,
      y=-1*log10(Cat_pval)), 
    fill = "firebrick",
    size = 1.5,
    shape = 23
  )  +
  geom_text_repel(
    data=dat_in_plot[which(dat_in_plot$category == "Obs" &
                             dat_in_plot$Median > dat_in_plot$top_9X &
                             dat_in_plot$Cat_pval < 1e-10),],
    aes(
      x=Median_pos,
      y=-1*log10(Cat_pval), 
      label = Symbol)
      )  +
  #geom_point(data=win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0][chr.x=="2L"],
  #           aes(x=(start/2 +end/2) , y=1*log10(rbinom.p)/10,
  #               color=(rnp.pr)), size=.95) +
  #scale_color_viridis(option="G", direction = -1) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~Analysis,
             ncol = 1) +
  ylab(expression( paste(italic(P[binomial])) ) ) +
  xlab("bp") +
  theme_bw() +
  ggtitle("Enrrichment of functional loci") +
  geom_hline(yintercept = -log10(10e-10)) + 
  geom_vline(xintercept = 2225744) + 
  geom_vline(xintercept = 13154180) -> 
  p_tresh_plot_genes

ggsave(p_tresh_plot_genes, 
       file = "p_tresh_plot_genes.pdf",
       width = 8,
       height = 8
       )

#### Enrichment
win.out.set = win.out[pr==.05][order(rnp.pr)][nSNPs>100][perm==0]
setDT(win.out.set)
setkey(win.out.set,  chr.x, start, end)

dat_in %<>%
  mutate(start=as.numeric(Median_pos),
           end=as.numeric(Median_pos),
           chr.x="2L")
setDT(dat_in)


foverlaps(dat_in[,c("chr.x","start","end","Analysis")], win.out.set[,c("start","end","chr.x","rbinom.p","invName")] , nomatch=NA) %>% 
  dplyr::select(Median_pos=i.start, Analysis, window_p=rbinom.p, invName) ->
  windows_p_correspondence
  
#save(dat_in, windows_p_correspondence, file = "checkpoint_window_vs_gene.Rdata" )
#load("./checkpoint_window_vs_gene.Rdata")
left_join(dat_in[which(dat_in$category == "Obs"),], windows_p_correspondence) %>% 
  .[complete.cases(.$invName),] ->
  dat_in_plus_win

group = unique(dat_in_plus_win$Analysis)

cors <- foreach(i=1:length(group), .combine = "rbind")%do%{
  dat_tmp = dat_in_plus_win[which(dat_in_plus_win$Analysis == group[i])]
  
  tmp_obj1 <- cor.test(~ window_p + Cat_pval, data = dat_tmp[which(dat_tmp$invName == "2Lt")])
  tmp_obj2 <- cor.test(~ window_p + Cat_pval+0.1, data = dat_tmp[which(dat_tmp$invName == "none")])
  
  rbind(
  data.frame(
    Analysis = group[i],
    cor = tmp_obj1$estimate,
    cor.p= tmp_obj1$p.value,
    type = "2lt"
  ),
  data.frame(
    Analysis = group[i],
    cor = tmp_obj2$estimate,
    cor.p= tmp_obj2$p.value,
    type = "none"))
  
}


dat_in_plus_win %>%
  mutate(beats_perm = case_when(Median >= top_9X ~ "yes",
                                Median < top_9X ~ "no") ) %>%
  ggplot(
    aes(
      x=-1*log10(Cat_pval),
      y=-1*log10(window_p),
      color = beats_perm
    )
  ) +
  geom_point() +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = -1*log10(0.05)) +
  geom_smooth(method = "lm", color = "black") +
  ylab("P value of window") +
  xlab("P value of gene") +
  facet_grid(invName~Analysis, scales = "free")->
  win_gene_plot

ggsave(win_gene_plot, file = "win_gene_plot.png",
       width = 7,
       height = 3)



  
  