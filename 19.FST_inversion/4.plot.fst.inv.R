###### Plot FST 
###### 
###### 
rm(list = ls())

library(tidyverse)
library(reshape2)
library(magrittr)
library(foreach)
library(viridis)
library(data.table)
library(ggridges)
library(broom)

###
files <- system("ls | grep 'Year_to_year_object.Rdata' ", intern = T)

###
fst.dat =
foreach(i=1:length(files), 
        .combine = "rbind")%do%{
  file.to.load = files[i]
  load(file.to.load) 
  Out_comp_vector_samepops
}

####
####
####
### -------> filter by effective coverage

inmeta="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata"
load(inmeta)

samps_EFFCOV %>% 
  filter(MeanEC > 28) %>% 
  .$sampleId -> samps.passECfilt

#### load demographic data
file.clusters = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/DEST_Sample_clusters.txt"
clust.asg = fread(file.clusters) %>%
  dplyr::select(samp1=sampleId, Continental_clusters)


####
####
####
####
fst.dat %>%
  left_join(clust.asg) %>%
  mutate(Continental_clusters = case_when(pop1 == "Charlottesville" ~ "Cville",
                                          TRUE ~ Continental_clusters)) %>%
  filter(samp1 %in% samps.passECfilt & samp2 %in% samps.passECfilt) ->
  fst.dat.EC

### Box plots
fst.dat.EC %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  mutate(month_diff = abs(month1-month2) ) %>% 
  ggplot(
    aes(x= as.factor(year_diff),
        y=FST,
        color = SNP.set)
  ) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  geom_boxplot(outlier.size =  0.5, size = 0.6) +
  facet_wrap(~Continental_clusters)->
  fst.inv.plot.boxplot

ggsave(fst.inv.plot.boxplot, file = "fst.inv.plot.boxplot.pdf", h = 2.5, w = 6)

### larger boxplot
fst.dat.EC %>%
  mutate(pop.set = case_when(pop1 %in% c("Akaa","Broggingen", "Odesa", "Charlottesville", "Munich") ~ pop1,
                             TRUE ~ Continental_clusters
                             )) %>%
  mutate(pop.set = factor(pop.set, levels = c("Charlottesville", "1.Europe_W", "Broggingen", "Munich", "Akaa", "3.Europe_E", "Odesa" ) ) ) %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  mutate(month_diff = abs(month1-month2) ) %>% 
  ggplot(
    aes(x= as.factor(year_diff),
        y=FST,
        color = SNP.set)
  ) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  geom_boxplot(outlier.size =  0.5, size = 0.6) +
  theme(legend.position = "bottom") +
  facet_grid(~pop.set)->
  fst.inv.plot.boxplot.poplevel

ggsave(fst.inv.plot.boxplot.poplevel, file = "fst.inv.plot.boxplot.poplevel.pdf", h = 3, w = 7)

#### Spatial FST
#### 
meta.for.merge = samps_EFFCOV
names(meta.for.merge)[1] = "samp1"

fst.dat.EC %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  filter(year_diff == 0) %>%  
  left_join(meta.for.merge, by = "samp1") %>%
  filter(Continental_clusters %in% c("1.Europe_W", "3.Europe_E"))  ->
  dat2015
  
dat2015 %>%
  group_by(pop1, lat, long, SNP.set, Continental_clusters) %>%
  summarise(mean.fst = mean(FST)) %>%
  ggplot(aes(
    x=lat,
    y=mean.fst,
    color = SNP.set 
    #shape = as.factor(year),
  )) +
  geom_jitter() +
  theme_bw() +
  ggtitle("In2Lt in Europe") +
  scale_color_viridis(discrete = TRUE, option = "D") +
  geom_smooth(method = "lm", se = F) +
  facet_grid(Continental_clusters ~ ., scales = "free_x") ->
  lat.reg.analysis

ggsave(lat.reg.analysis, file = "lat.reg.analysis.pdf", h = 3, w = 5)

dat2015 %>%
  group_by(pop1, lat, long, SNP.set, Continental_clusters) %>%
  summarise(mean.fst = mean(FST)) %>%
  ggplot(aes(
    x=long,
    y=mean.fst,
    color = SNP.set 
    #shape = as.factor(year),
  )) +
  geom_jitter() +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  ggtitle("In2Lt in Europe") +
  geom_smooth(method = "lm", se = F) +
  facet_grid(Continental_clusters ~ ., scales = "free_x") ->
  lon.reg.analysis

ggsave(lon.reg.analysis, file = "lon.reg.analysis.pdf", h = 3, w = 5)

cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "3.Europe_E" & SNP.set == "glm.snps"  ) )
cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "3.Europe_E" & SNP.set == "macthed.controls.Inv"  ) )
cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "3.Europe_E" & SNP.set == "macthed.controls.noInv"  ) )

cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "1.Europe_W" & SNP.set == "glm.snps"  ) )
cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "1.Europe_W" & SNP.set == "macthed.controls.Inv"  ) )
cor.test( ~ FST + lat , data = filter(dat2015, Continental_clusters == "1.Europe_W" & SNP.set == "macthed.controls.noInv"  ) )

###
fst.dat.EC %>%
  filter(pop1 %in% c("Charlottesville", "Odesa")  ) %>%
  mutate(year_diff = abs(year1-year2) ) %>% 
  mutate(month_diff = abs(month1-month2) ) %>% 
  group_by(SNP.set, year_diff,month_diff, pop1#, year1, year2
           ) %>%
  filter(year_diff %in% 0:1) %>% 
  summarise(median.FST = quantile(FST, 0.5),
            mean.FST = mean(FST),
            ) %>%
  ggplot(
    aes(
     x= month_diff,
     y= mean.FST,
    color = SNP.set
  )) +
  #geom_density_ridges(scale = 4, adjust = 100) +
  #geom_smooth(se = F) +
  #geom_boxplot() +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, option = "D") +
  facet_grid(as.factor(year_diff)~pop1) ->
  #facet_grid(year1~ year2) ->
  fst.inv.plot.mobth

ggsave(fst.inv.plot.mobth, file = "fst.inv.plot.mobth.pdf",h = 2.5, w = 6)

  
  ####
  ####
  ####
  ####
### add weather data
  ######
  load("/scratch/yey2sn/Overwintering_ms/16.Haplotypes/nasa_power.weather.mod.Rdata")
  names(weather.ave)[1] = "sampleId"
  
  sets <- data.table(mod=c(1:11),
                     start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                     end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
  
  weather.ave %>%
    filter(mod == 2) %>%
    dplyr::select(sampleId, temp.max) -> tmax
  
  weather.ave %>%
    filter(mod == 9) %>%
    dplyr::select(sampleId, temp.ave) -> tave
  
  weather.ave %>%
    filter(mod == 8) %>%
    dplyr::select(sampleId, humidity.ave) -> have
  
  cbind(tmax,  temp.ave=tave[,-1], humid.ave=have[,-1]) %>%
    as.data.frame() %>%
    group_by(sampleId) %>%
    slice_head(n=1) -> obs_for_cors

  obs_for_cors %>% dplyr::select(samp1 = sampleId, temp1 = temp.max) -> t1
  obs_for_cors %>% dplyr::select(samp2 = sampleId, temp2 = temp.max) -> t2
  
  obs_for_cors %>% dplyr::select(samp1 = sampleId, temp1 = temp.ave.temp.ave) -> ta1
  obs_for_cors %>% dplyr::select(samp2 = sampleId, temp2 = temp.ave.temp.ave) -> ta2
  
  obs_for_cors %>% dplyr::select(samp1 = sampleId, temp1 = humid.ave.humidity.ave) -> h1
  obs_for_cors %>% dplyr::select(samp2 = sampleId, temp2 = humid.ave.humidity.ave) -> h2
  
  
  ######
  fst.dat.EC %>%
    filter(Continental_clusters == "Cville") %>%
    left_join(t1) %>% 
    left_join(t2) %>% 
    mutate(temp_diff = abs(temp1-temp2)) %>%
    mutate(year_diff = abs(year1-year2) ) %>% 
    mutate(month_diff = abs(month1-month2) ) %>% 
    group_by(SNP.set,month_diff, year_diff,temp_diff, #, year1, year2
    ) %>%
    filter(year_diff %in% 0:1) %>%
    summarise(median.FST = quantile(FST, 0.5),
              mean.FST = mean(FST),
    )  %>% mutate(clust = "Cville") -> cville.temp
   
  fst.dat.EC %>%
    filter(pop1 == "Odesa") %>%
    left_join(ta1) %>% 
    left_join(ta2) %>% 
    mutate(temp_diff = abs(temp1-temp2)) %>%
    mutate(year_diff = abs(year1-year2) ) %>% 
    mutate(month_diff = abs(month1-month2) ) %>% 
    group_by(SNP.set,month_diff, year_diff,temp_diff, #, year1, year2
    ) %>%
    filter(year_diff %in% 0:1) %>%
    summarise(median.FST = quantile(FST, 0.5),
              mean.FST = mean(FST),
    )  %>% mutate(clust = "Odesa") -> odesa.temp
  
   
  fst.dat.EC %>%
    filter(Continental_clusters == "1.Europe_W") %>%
    left_join(h1) %>% 
    left_join(h2) %>% 
    mutate(temp_diff = abs(temp1-temp2)) %>%
    mutate(year_diff = abs(year1-year2) ) %>% 
    mutate(month_diff = abs(month1-month2) ) %>% 
    group_by(SNP.set,month_diff, year_diff,temp_diff, #, year1, year2
    ) %>%
    filter(year_diff %in% 0:1) %>%
    summarise(median.FST = quantile(FST, 0.5),
              mean.FST = mean(FST),
    ) %>% mutate(clust = "EUW") -> euw.hu
  
  fst.dat.EC %>%
    filter(Continental_clusters == "3.Europe_E") %>%
    left_join(ta1) %>% 
    left_join(ta2) %>% 
    mutate(temp_diff = abs(temp1-temp2)) %>%
    mutate(year_diff = abs(year1-year2) ) %>% 
    mutate(month_diff = abs(month1-month2) ) %>% 
    group_by(SNP.set,month_diff, year_diff,temp_diff, #, year1, year2
    ) %>%
    filter(year_diff %in% 0:1) %>%
    summarise(median.FST = quantile(FST, 0.5),
              mean.FST = mean(FST),
    ) %>% mutate(clust = "EUE") -> eue.temp
  
  ## test
  
  rbind(cville.temp, 
        euw.temp, 
        eue.hu, 
        odesa.temp
  ) -> dat.for.model
  
  dat.for.model %>% nest(data = -c(SNP.set,year_diff,clust) ) %>% 
    mutate(model = map(data, ~lm(mean.FST~temp_diff, data = .)), tidied = map(model, tidy)) %>%
    unnest(tidied) %>%
    filter(term == "temp_diff") %>%
    dplyr::select(SNP.set, clust, year_diff, estimate, p.value, std.error) %>% 
    mutate(estimate = as.numeric(format(estimate, format = "e", digits = 2))) -> linear.models.fst
  
  linear.models.fst %>%
    ggplot(
      aes(
        x=clust,
        y=estimate,
        ymin=estimate-std.error,
        ymax=estimate+std.error,
        color= SNP.set
      )
    ) +
    geom_hline(yintercept = 0) +
    geom_errorbar(position=position_dodge(width=0.5), width = 0.1) +
    geom_point(position=position_dodge(width=0.5), size = 2.0) +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "D") +
    facet_grid(year_diff~.)->
    betas.fst
  
  ggsave(betas.fst, file = "betas.fst.pdf", h = 2.5, w = 4.5)
  
  ## plot
  
  rbind(cville.temp, 
        euw.temp, 
        eue.hu, 
        odesa.temp
        ) %>%
    ggplot(aes(
      x= temp_diff,
      y= mean.FST,
      color = SNP.set
    )) +
    #geom_boxplot() +
    #geom_line() +
    geom_point(alpha = 0.3, 
               size = 1.5, shape = 21, color = "black", aes(fill = SNP.set)) +
    geom_smooth(se = F, 
                #span = 10,
                method = "lm"
                ) +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    facet_grid(as.factor(year_diff) ~ clust, #scales = "free_y"
               ) ->
    #facet_grid(year1~ year2) ->
    fst.inv.plot.temp
  
  ggsave(fst.inv.plot.temp, file = "fst.inv.plot.temp.pdf", h = 2.5, w = 8)
  
  ####
   
  #### GEO ANALYSIS
  ####   #### GEO ANALYSIS
  #### GEO ANALYSIS
  #### GEO ANALYSIS
  #### GEO ANALYSIS
  #### GEO ANALYSIS
  #### GEO ANALYSIS
  #### GEO ANALYSIS
  
  ###
  files.geo <- system("ls | grep 'geo' ", intern = T)
  
  ###
  fst.dat.geo =
    foreach(i=1:length(files.geo), 
            .combine = "rbind")%do%{
              file.to.load = files.geo[i]
              load(file.to.load) 
              Out_comp_vector_samepops
            }
  
  
  inmeta="/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/DEST_EC_metadata.Rdata"
  load(inmeta)
  
  samps_EFFCOV %>% 
    filter(MeanEC > 28) %>% 
    .$sampleId -> samps.passECfilt
  
  #### load demographic data
  file.clusters = "/scratch/yey2sn/Overwintering_ms/19.inv.fst/DEST_Sample_clusters.txt"
  clust.asg = fread(file.clusters) %>%
    dplyr::select(samp1=sampleId, Continental_clusters)
  
  fst.dat.geo %>%
    left_join(clust.asg) %>%
    mutate(Continental_clusters = case_when(pop1 == "Charlottesville" ~ "Cville",
                                            TRUE ~ Continental_clusters)) %>%
    filter(samp1 %in% samps.passECfilt & samp2 %in% samps.passECfilt) ->
    fst.dat.geo.EC
  
  fst.dat.geo.EC %>%
    ggplot(aes(
      x=log10(geo.dist/1000),
      y=FST,
      color = SNP.set
    )) +
    geom_point() +
    theme_bw() +
    geom_smooth(method = "lm", se = F) +
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    facet_grid(~Continental_clusters) ->
    geo.fst
  
  ggsave(geo.fst, file = "geo.fst.pdf", h = 2.5, w = 8)
 
  fst.dat.geo.EC %>%
    ggplot(aes(
      x=Continental_clusters,
      y=FST,
      color = SNP.set
    )) +
    geom_boxplot() +
    theme_bw() +
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE, option = "D")  ->
    geo.fst.box
  
  ggsave(geo.fst.box, file = "geo.fst.box.pdf", h = 2.5, w = 6)
  
  
  