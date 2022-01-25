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

## set some generalities
p_tresh=0.05

##
##


######## ---> Year MODEL
######## 
### This line finsd the number of LRT p values of the model for real data
## Load data:
i=0
model="year_factor"

load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))

glm.out %>%
  filter(mod == model) %>%
  filter(p.lrt < p_tresh) %>% 
  dim %>% .[1] -> out_alph1

glm.out %>%
  filter(mod == model) %>% dim %>% .[1] -> out_all

out_alph1/out_all -> perc_real


### This line finsd the number of LRT p values of the model for the permutated data

percent_list= list()
for(i in 1:100){
  
  if(i == 99){break}
  
load(paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", i, ".Rdata", sep = ""))

glm.out %>%
filter(mod == model) %>%
filter(p.lrt < p_tresh) %>% 
  dim %>% .[1] -> out_alph1

glm.out %>%
filter(mod == model) %>% dim %>% .[1] -> out_all

out_alph1/out_all -> perc_perm
percent_list[[i]] = perc_perm
} # close i

print(paste("the mean % of time SNPS in perms is",
mean(unlist(percent_list)*100)))
print(paste("the SD % of time SNPS in perms is",
            sd(unlist(percent_list)*100)))


############ ########## ############
############ Summarize data from the time + tempearture model
############ 

#### Files were processed with the script: "Run_YT_model_pvalue_extractor.R"

pops = c(
  "VA_ch",
  "UA_Ode",
  "TR_Yes",
  "PA_li",
  "FI_Aka",
  "DE_Mun",
  "DE_Bro"
)

collect_counts = list()
collect_enrrich = list()
for(i in 1:length(pops)){
  print(pops[i])
  
  load(paste("/scratch/yey2sn/Overwintering_ms/4.GML_plots/",
             pops[i],
             ".enrrich_count.Rdata",
             sep = ""))
  
  collect_counts[[i]] = count_df
  collect_enrrich[[i]] = enrrich_df
  
}

### Part 1 -- counts
counts_df = do.call(rbind, collect_counts)
counts_df %<>%
  mutate(proportion = Tresh/N)

counts_df$chr[grep("genome", counts_df$chr)] = "G"

  ggplot() +
    geom_violin(data=counts_df[which(counts_df$category == "Per"),],
                aes(x=chr,
                    y=proportion),
                alpha = 0.7,
                fill = "springgreen") +
    geom_point(data=counts_df[which(counts_df$category == "Obs"),],
               aes(x=chr,
                   y=proportion),
               size = 4,
               shape = 23,
               color = "black",
               fill = "firebrick") +
    theme_bw() +
      facet_grid(~pop)  -> 
      p_tresh_plot

ggsave(p_tresh_plot, 
       file = "p_tresh_plot.pdf",
       width = 8,
       height = 2)


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


###
enrrich_df %>% filter(analysis_type == "enrr_Inv") -> 
  enrr_Inv

### prepare table
enrr_Inv %>%
  group_by(chr, inversion_pos, category, pop) %>%
  summarise(Nobs = median(N)) %>%
  dcast(chr+pop~inversion_pos+category, 
        fill = 0,
        value.var = "Nobs") ->
  enrr_Inv_cast

inv_enrrichment = list()
for(i in 1:dim(enrr_Inv_cast)[1]){
  tmp_matrix <-
  matrix(rep(NA, 4),
         nrow = 2,
         dimnames = list( c("Obs", "Perm"),
                          c("Feature", "Not_Feature") ) )

tmp_matrix[1,1] =  enrr_Inv_cast[i,1+2]
tmp_matrix[1,2] =  enrr_Inv_cast[i,3+2]
tmp_matrix[2,1] = enrr_Inv_cast[i,2+2]
tmp_matrix[2,2] = enrr_Inv_cast[i,4+2]

fet_tmp <- fisher.test(tmp_matrix)

tmp_inv = data.frame(
  type = "inv",
  chr = enrr_Inv_cast[i, "chr"],
  pop = enrr_Inv_cast[i, "pop"],
  OR = fet_tmp$estimate,
  OR_low = fet_tmp$conf.int[1],
  OR_high = fet_tmp$conf.int[2],
  OR_p = fet_tmp$p.value
)

inv_enrrichment[[i]] = tmp_inv
}

inv_enrrichment_df = do.call(rbind, inv_enrrichment)
inv_enrrichment_df$OR[inv_enrrichment_df$chr == "2R"] = 1
inv_enrrichment_df$OR_low[inv_enrrichment_df$chr == "2R"] = 1
inv_enrrichment_df$OR_high[inv_enrrichment_df$chr == "2R"] = 1

inv_enrrichment_df %>%
  ggplot(aes(
    x=chr,
    y=log2(OR),
    ymin =log2(OR_low),
    ymax =log2(OR_high),
    fill = -log10(OR_p)
  )) +
  geom_hline(yintercept = log2(1), linetype = "dashed") + 
  geom_errorbar(width = 0.01) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradient2(mid = -log10(0.01)) +
  ggtitle("Enrrichment test Inversion vs. No Inversion") + 
  facet_grid(~pop) ->
  p_tresh_enrr_Inv

ggsave(p_tresh_enrr_Inv, 
       file = "p_tresh_enrr_Inv.pdf",
       width = 9,
       height = 3)

#######

enrrich_df %>% filter(analysis_type == "enrr_consq") -> 
  enrr_Consq

### prepare table
enrr_Consq %>%
  group_by(chr, inversion_pos, category, pop, Consequence) %>%
  summarise(Nobs = median(N)) %>%
  dcast(chr+pop+inversion_pos+Consequence~category, 
        fill = 0,
        value.var = "Nobs") %>%
  filter(Consequence %in% c("3_prime",
                            "5_prime",
                            "intergenic_region",
                            "intron_variant",
                            "missense_variant",
                            "Splice_var",
                            "synonymous_variant")) ->
  enrr_Consq_cast

enrr_Consq_cast %>% 
  group_by(chr, pop) %>%
  summarise(N_obs = sum(Obs),
            N_per = sum(Per)) ->
  total_obs_table

left_join(enrr_Consq_cast, 
          total_obs_table) %>%
  mutate(Obs_Rest = N_obs-Obs,
         Per_Rest = N_per-Per) ->
  enrrichment_conseq_in_data


consq_enrrichment = list()
for(i in 1:dim(enrrichment_conseq_in_data)[1]){
  
  tmp_matrix <-
    matrix(rep(NA, 4),
           nrow = 2,
           dimnames = list( c("Obs", "Perm"),
                            c("Feature", "Not_Feature") ) )
  
  tmp_matrix[1,1] =  enrrichment_conseq_in_data[i, "Obs"]
  tmp_matrix[1,2] =  enrrichment_conseq_in_data[i, "Obs_Rest"]
  tmp_matrix[2,1] =  enrrichment_conseq_in_data[i, "Per"]
  tmp_matrix[2,2] =  enrrichment_conseq_in_data[i, "Per_Rest"]
  
  fet_tmp <- fisher.test(tmp_matrix)
  
  tmp_cons = data.frame(
    type = "conseq",
    Conseq = enrrichment_conseq_in_data[i, "Consequence"],
    chr = enrrichment_conseq_in_data[i, "chr"],
    pop = enrrichment_conseq_in_data[i, "pop"],
    inv = enrrichment_conseq_in_data[i, "inversion_pos"],
    OR = fet_tmp$estimate,
    OR_low = fet_tmp$conf.int[1],
    OR_high = fet_tmp$conf.int[2],
    OR_p = fet_tmp$p.value
  )
  
  consq_enrrichment[[i]] = tmp_cons
}

consq_enrrichment_df = do.call(rbind, consq_enrrichment)

consq_enrrichment_df %>%
  mutate(signif = case_when(OR_p <= 0.05 ~ "S",
                            OR_p > 0.05 ~ " NSig" )) %>%
  ggplot(aes(
    x=chr,
    y=log2(OR),
    ymin =log2(OR_low),
    ymax =log2(OR_high),
    color = signif,
    shape = inv,
    group = inv
  )) +
  geom_hline(yintercept = log2(1), linetype = "dashed") + 
  geom_errorbar(width = 0.01, position=position_dodge(width=0.5) ) +
  geom_point(size = 2, position=position_dodge(width=0.5) ) +
  #scale_fill_gradient2(mid = -log10(0.01)) +
  ylim(-1,1) + 
  ggtitle("Enrrichment test Consequence") + 
  facet_grid(Conseq~pop) ->
  p_tresh_enrr_Con

ggsave(p_tresh_enrr_Con, 
       file = "p_tresh_enrr_Con.pdf",
       width = 8,
       height = 9)

