## summarizing SNPs

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

###

library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(gmodels)

############ ########## ############
############ Summarize data from the time + tempearture model
############ 

#### Files were processed with the script: "Run_YT_model_pvalue_extractor.R"

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

counts_df$analysis_subtype[grep("genome", counts_df$analysis_subtype)] = "G"

counts_df %<>%
  mutate(fill_for_plot = paste(inversion_pos,CRM_pos, sep = "")) 
counts_df$fill_for_plot = gsub("NA", "", counts_df$fill_for_plot )
counts_df$fill_for_plot = gsub("^$", "All", counts_df$fill_for_plot )
counts_df$fill_for_plot %>% table

### Part 0 -- Exploring the distribution of counts as a fucntion of p.values
Cohens_h_f = function(p1, p2){
  phi1 = 2*asin(sqrt(p1))
  phi2 = 2*asin(sqrt(p2))
  return(phi1-phi2)
}

counts_df %>%
  group_by(category,p_tresh,pop, fill_for_plot, chr) %>%
  summarise(Med_p = median(proportion)) %>% 
  #dcast(chr+fill_for_plot+p_tresh+pop~category) %>% 
  #mutate(Cohens_h = Cohens_h_f(Obs,Per)) %>% 
  ggplot(aes(
         x= as.numeric(p_tresh),
         #y= (Cohens_h),
         y=Med_p,
         color = chr,
         linetype = category
        )) +
  geom_line() +
  #coord_trans(x = "log10", y = "log10") +
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text = element_text(size = 7, angle=45)) +
  ylab( expression(log[10]("%SNPs" < alpha)) ) +
  xlab(expression(log[10](alpha))) +
  facet_grid(pop~fill_for_plot) -> p.value.by.count.plot

ggsave(p.value.by.count.plot, 
       file = "p.value.by.count.plot.pdf",
       width = 12,
       height = 8)

### Examine distribution of P-values

### Part 1 -- Digest the count distributions

dat_in = counts_df %>%
  filter(#pop == "VA_ch",
         p_tresh == 0.05,
         fill_for_plot %in% c("All", "inv", "non.inv")) 
  
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
    theme_bw() +
    theme(legend.position = "top")-> 
      p_tresh_plot

ggsave(p_tresh_plot, 
       file = "p_tresh_plot.pdf",
       width = 9,
       height = 2.5)


### Statistics of the dataset

#pops
#counts_df$chr[grep("genome", counts_df$chr)] = "G"

chrs = c("G", "2L", "2R", "3L", "3R")

for(i in 1:length(pops) ) {
  for(j in 1:length(chrs) ){
    
    Ho <- counts_df$proportion[which(counts_df$pop == pops[i] & counts_df$chr ==chrs[j] & counts_df$category == "Obs")] 
    DatProp <- counts_df$proportion[which(counts_df$pop == pops[i] & counts_df$chr ==chrs[j] & counts_df$category == "Per")] 
    
    ecdf(DatProp) -> tmp_ecdf
    
    1-tmp_ecdf(Ho)
    
    data.frame()
    
  } # close j
} #close i

Ho <- counts_df$proportion[which(counts_df$pop == "VA_ch" & counts_df$chr =="G" & counts_df$category == "Obs")] 
DatProp <- counts_df$proportion[which(counts_df$pop == "VA_ch" & counts_df$chr =="G" & counts_df$category == "Per")] 
quantile(DatProp, Ho)


### Part 2 -- enrrichment
enrrich_df = do.call(rbind, collect_enrrich)

## add the CRMS
##counts_df %>% 
##  .[complete.cases(.),] %>% 
##  filter(CRM_pos == "CRM") %>% 
##  as.data.frame() %>%
##  dplyr::select(chr,
##         Consequence = CRM_pos,
##         inversion_pos,
##         N = Tresh,
##         p_tresh,
##         category,
##         pop) %>% 
##  mutate(analysis_type = "enrr_consq") -> CRMs_counts_enrr
##
#######

enrrich_df %>% filter(analysis_type == "enrr_consq") %>%
filter(p_tresh >= 0.01) -> enrr_Consq

##rbind(CRMs_counts_enrr, enrr_Consq) -> enrr_Consq

### prepare table
enrr_Consq %>%
  group_by(chr, category, pop, Consequence, p_tresh) %>%
  filter(pop == "VA_ch") %>% 
  filter(category == "Obs") -> obs_only

names(obs_only)[3] = "Obs_N"

enrr_Consq %>%
  group_by(chr, inversion_pos, category, pop, Consequence, p_tresh) %>%
  filter(pop == "VA_ch") %>% 
  filter(category == "Per") -> Perm_only

names(Perm_only)[3] = "Per_N"


left_join(Perm_only[-5], obs_only[-5]) %>%
  #summarise(Nobs = median(N)) %>%
  #dcast(chr+p_tresh+pop+inversion_pos+Consequence~category, 
        #fill = 0,
        #value.var = "Nobs") %>%
  filter(!Consequence %in% c("NC_exon",
                             "start",
                             "stop")) ->
  enrr_Consq_cast

enrr_Consq_cast %>% 
  group_by(chr, pop, p_tresh, inversion_pos) %>% 
  summarise(Obs_tot = sum(Obs_N),
            Per_tot = sum(Per_N)) ->
  total_obs_table

left_join(enrr_Consq_cast, 
          total_obs_table) %>%
  mutate(Obs_Rest = Obs_tot-Obs_N,
         Per_Rest = Per_tot-Per_N) %>%
  as.data.frame() %>%
  .[complete.cases(.),]->
  enrrichment_conseq_in_data


##### Sanity cheack

dat_in2 = enrrichment_conseq_in_data %>%
  mutate(Tidy_Name = case_when(
    Consequence == "3_prime" ~ "3p",
    Consequence == "5_prime" ~ "5p",
    Consequence == "intergenic_region|CRM" ~ "Inter+CRM",
    Consequence == "intergenic_region|non.CRM" ~ "Inter",
    Consequence == "intron_variant|CRM" ~ "Intron+CRM",
    Consequence == "intron_variant|non.CRM" ~ "Intron",
    Consequence == "missense_variant" ~ "Nonsyn",
    Consequence == "synonymous_variant" ~ "Syn"
  )) %>% 
  .[complete.cases(.$Tidy_Name),] %>%
  filter(pop == "VA_ch",
         analysis_type == "enrr_consq",
         p_tresh == 0.05) %>%
  mutate(Obs_prop = Obs_N/Obs_tot,
         Per_prop = Per_N/Per_tot)  %>%
  dplyr::select(chr, inversion_pos, Tidy_Name,   Obs_prop ,  Per_prop) %>% 
  melt(id = c("chr", "inversion_pos", "Tidy_Name")) %>% 
  separate(variable, into = c("category","metric")) 
  
dat_in2 %>% 
  filter(category == "Per") -> dat_per

dat_in2 %>% 
  filter(category == "Obs") %>% 
  group_by(chr, inversion_pos, Tidy_Name, category, metric) %>% 
  summarise(value = median(value)) -> dat_obs

dat_in2 = rbind(dat_obs, dat_per)

ggplot() +
  geom_violin(data=dat_in2[which(dat_in2$category == "Per"),],
              aes(x=chr,
                  y=value,
                  fill = inversion_pos),
              alpha = 0.7
  ) +
  geom_point(data=dat_in2[which(dat_in2$category == "Obs"),],
             aes(x=chr,
                 y=value,
                 fill = inversion_pos),
             size = 2.3,
             shape = 23,
             color = "black",
             position = position_dodge2(w = 0.95)) +
  facet_grid(~Tidy_Name, scales = "free_x", space = "free") +
  ylab(expression( paste("%", italic(P[LRT]<0.05)) ) ) +
  theme_bw() +
  theme(legend.position = "top")-> 
  p_tresh_plot

ggsave(p_tresh_plot, 
       file = "p_tresh_plot.pdf",
       width = 9,
       height = 2.5)



#####

consq_enrrichment = list()
for(i in 1:dim(enrrichment_conseq_in_data)[1]){
  
  print(i/dim(enrrichment_conseq_in_data)[1]*100)
  
  tmp_matrix <-
    matrix(rep(NA, 4),
           nrow = 2,
           dimnames = list( c("Obs", "Perm"),
                            c("Feature", "Not_Feature") ) )
  
  tmp_matrix[1,1] =  enrrichment_conseq_in_data[i, "Obs_N"]
  tmp_matrix[1,2] =  enrrichment_conseq_in_data[i, "Obs_Rest"]
  tmp_matrix[2,1] =  enrrichment_conseq_in_data[i, "Per_N"]
  tmp_matrix[2,2] =  enrrichment_conseq_in_data[i, "Per_Rest"]
  
  fet_tmp <- fisher.test(tmp_matrix)
  
  tmp_cons = data.frame(
    type = "conseq",
    Conseq = enrrichment_conseq_in_data[i, "Consequence"],
    chr = enrrichment_conseq_in_data[i, "chr"],
    pop = enrrichment_conseq_in_data[i, "pop"],
    inv = enrrichment_conseq_in_data[i, "inversion_pos"],
    p_tresh = enrrichment_conseq_in_data[i, "p_tresh"],
    OR = fet_tmp$estimate,
    OR_low = fet_tmp$conf.int[1],
    OR_high = fet_tmp$conf.int[2],
    OR_p = fet_tmp$p.value
  )
  
  consq_enrrichment[[i]] = tmp_cons
}

consq_enrrichment_df = do.call(rbind, consq_enrrichment)

save(consq_enrrichment_df, file = "./consq_enrrichment_df.Rdata")

load("./consq_enrrichment_df.Rdata")

##consq_enrrichment_df %>%
##  mutate(signif = case_when(OR_p <= 0.05 ~ "S",
##                            OR_p > 0.05 ~ " NSig" )) %>%
##  filter(
##         p_tresh >= 0.01) %>%
##  ggplot(aes(
##    x=as.numeric(p_tresh),
##    y=log2(OR),
##    #ymin =log2(OR_low),
##    #ymax =log2(OR_high),
##    color = Conseq,
##    #alpha = -log10(OR_p)
##    #shape = inv,
##    #group = inv
##  )) +
##  geom_hline(yintercept = log2(1), linetype = "dashed") + 
##  geom_point() +
##  #geom_errorbar(width = 0.01, position=position_dodge(width=0.5) ) +
##  #geom_point(size = 2, position=position_dodge(width=0.5) ) +
##  #scale_fill_gradient2(mid = -log10(0.01)) +
##  ylim(-1,1) + 
##  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
##  labels = trans_format("log10", math_format(10^.x)),) +
##  ggtitle("Enrrichment test Consequence, ALL POPS") + 
##  theme(axis.text = element_text(size = 7, angle=45)) +
##  ylab( expression(log[2](OR)) ) +
##  xlab(expression(log[10](alpha))) +
##  facet_grid(chr+inv~pop) ->
##  p_tresh_enrr_Con
##
##ggsave(p_tresh_enrr_Con, 
##       file = "p_tresh_enrr_Con.pdf",
##       width = 12,
##       height = 9)

### VA plot
consq_enrrichment_df %>%
  mutate(Tidy_Name = case_when(
    Conseq == "3_prime" ~ "3p",
    Conseq == "5_prime" ~ "5p",
    Conseq == "intergenic_region|CRM" ~ "Inter+CRM",
    Conseq == "intergenic_region|non.CRM" ~ "Inter",
    Conseq == "intron_variant|CRM" ~ "Intron+CRM",
    Conseq == "intron_variant|non.CRM" ~ "Intron",
    Conseq == "missense_variant" ~ "Nonsyn",
    Conseq == "synonymous_variant" ~ "Syn"
  )) %>% 
  #filter(OR_p <= 0.05) %>%
  #mutate(OR_corrected = case_when(p.adjust(OR_p) <= 0.05 ~ OR,
  #p.adjust(OR_p) > 0.05 ~ 1 )) %>% 
  #separate(Conseq, into = c("Conseq","CRM"), sep = "\\|") %>% 
  filter(pop == "VA_ch",
         p_tresh >= 0.01,
         !Conseq == "Splice_var") %>%
  group_by(Conseq, chr, pop, inv, p_tresh, Tidy_Name) %>%
  summarise(OR_mean = mean(OR),
            OR_sd = sd(OR),
            OR_l = ci(OR, confidence=0.99)[2],
            OR_h = ci(OR, confidence=0.99)[3],
            OR_lowQ = quantile(OR, 0.025),
            OR_highQ = quantile(OR, 0.975),
            ) %>% 
  ##consq_enrrichment_df_sep$CRM[is.na(consq_enrrichment_df_sep$CRM)] = "non.CRM"
  #consq_enrrichment_df_sep %>%
  ggplot(aes(
    x=as.numeric(p_tresh),
    y=log2(OR_mean),
    ymin =log2(OR_lowQ),
    ymax =log2(OR_highQ),
    fill = inv,
    color = inv
    #linetype = CRM
    #alpha = -log10(OR_p)
    #shape = inv,
    #group = inv
  )) +
  geom_hline(yintercept = log2(1), linetype = "dashed") + 
  #geom_jitter(size = 0.5, alpha = 0.1) +
  #geom_smooth() +
  #geom_boxplot(outlier.shape = NA) +
  geom_ribbon(alpha = 0.1) +
  geom_line() +
  #geom_errorbar(width = 0.01, position=position_dodge(width=0.1) ) +
  #geom_point(size = 1, position=position_dodge(width=0.1) ) +
  #scale_fill_gradient2(mid = -log10(0.01)) +
  #ylim(-0.8,0.8) +
  #coord_trans(x="log10") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),) +
  ggtitle("Enrrichment test of Charlottesville, VA") + 
  theme(axis.text = element_text(size = 7, angle=45)) +
  ylab( expression(log[2](OR)) ) +
  xlab(expression(log[10](alpha))) +
  facet_grid(chr~Tidy_Name) ->
  p_tresh_enrr_Con_VA

ggsave(p_tresh_enrr_Con_VA, 
       file = "p_tresh_enrr_Con_VA.pdf",
       width = 12,
       height = 6)

####

ptresh = c(0.01,  0.05)

for(i in 1:2){
  consq_enrrichment_df %>%
  filter(p_tresh == ptresh[i],
         !Conseq == "Splice_var",
         OR_p < 0.05) %>%
  group_by(Conseq, chr, inv, p_tresh)  %>%
  summarise(OR_mean = ci(OR)[1],
            OR_l = ci(OR)[2],
            OR_h = ci(OR)[3]) %>%
  ggplot(aes(
    x=chr,
    y=log2(OR_mean),
    ymin =log2(OR_l),
    ymax =log2(OR_h),
    fill = inv,
    color = inv,
    #linetype = CRM,
    #shape = signif
    #alpha = -log10(OR_p)
    #shape = inv,
    #group = inv
  )) +
  geom_hline(yintercept = log2(1), linetype = "dashed") + 
  #geom_ribbon(alpha = 0.1) +
  #geom_line() +
  #geom_boxplot() +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5) ) +
  geom_point(size = 2, position=position_dodge(width=0.5) ) +
  #scale_fill_gradient2(mid = -log10(0.01)) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)),) +
  ggtitle(paste("alpha = ", ptresh[i],"Enrrichment test Consequence, VA") ) + 
  theme(axis.text = element_text(size = 7, angle=45)) +
  ylab( expression(log[2](OR)) ) +
  xlab(expression(log[10](alpha))) +
  facet_grid(~Conseq) ->
  p_tresh_5per_enrr_Con_VA

ggsave(p_tresh_5per_enrr_Con_VA, 
       file = paste(ptresh[i], "p_tres_enrr_Con_VA.pdf", sep = "."),
       width = 12,
       height = 4)
}

