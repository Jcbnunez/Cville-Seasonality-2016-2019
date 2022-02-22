## Haplotype analysis
## 
## 
## 
## Part 1 -- Identify inversion status

library(tidyverse)
library(magrittr)
library(data.table)

## Load SVM predictions
## uses ../Cville-Seasonality-2016-2019/5.Finding_inv2Lt_markers/9.train_predictive_model.r

load("/project/berglandlab/Dmel_genomic_resources/Inversions/2LT_inversion_markers/SVM_2ltpred_model_and_Files.Rdata")

CM_predictions %>%
  ggplot(aes(SVM)) +
  geom_vline(xintercept = 0.90) +
  geom_vline(xintercept = 0.10) +
  geom_histogram() -> svm_hist

ggsave(svm_hist, file = "svm_hist.pdf")

CM_predictions %>%
  separate(line_name, into = c("set", "id", "date"), sep = "_") %>%
  mutate(date_parsed = case_when( date %in% c("CMspring") ~ as.Date("06/01/2018", format = "%m/%d/%Y"),
                                  date %in% c("CMfall") ~ as.Date("09/01/2017", format = "%m/%d/%Y"),
                                  !(date %in% c("CMspring", "CMfall")) ~ as.Date(paste(date, 2016, sep = ""), format = "%m%d%Y"))) %>%
  mutate(karyot = case_when(SVM >= 0.90 ~ "Inv/Inv",
                            SVM <= 0.10 ~ "Std/Std",
                            SVM > 0.10 | SVM < 0.90 ~ "Std/Inv" ) ) %>%
  mutate(inv_freq = case_when(karyot == "Inv/Inv" ~ 0,
                              karyot == "Std/Std" ~ 2,
                              karyot == "Std/Inv" ~ 1) ) %>%
  mutate(geno = case_when(karyot == "Inv/Inv" ~ "Inv",
                              karyot == "Std/Std" ~ "Std",
                              karyot == "Std/Inv" ~ "Het") ) %>%
  mutate(Group = case_when(set == "CM" ~ "Alys",
                           set != "CM" ~ "Pris",) ) ->
  CM_predictions_parsed


CM_predictions_parsed %>% 
  group_by(date_parsed) %>%
  summarise(chr_tod = n()*2) -> tots


CM_predictions_parsed %>% 
  group_by(date_parsed,  Group) %>%
  summarise(gen_count = sum(inv_freq)) -> counts

left_join(tots, counts ) %>%
  mutate(prop = NA,
         prop_lo = NA,
         prop_h = NA) -> prop_table

for(i in 1:dim(prop_table)[1]){
  
  prop.test(prop_table$gen_count[i], 
            n = prop_table$chr_tod[i]) -> tmp

  prop_table$prop[i] = tmp$estimate
  prop_table$prop_lo[i] = tmp$conf.int[1]
  prop_table$prop_h[i] = tmp$conf.int[2]
  
}

prop_table %>%
  ggplot(aes(x=date_parsed, y = prop, 
             #ymin =prop_lo, 
             ##ymax = prop_h
             )) +
  #geom_errorbar() +
  geom_line() +
  geom_point() +
  ggtitle("Frequency of the STD karyotypes in Alys/Pris Datasets") +
  facet_wrap(~Group, scales ="free", ncol = 1)-> karyot_time

ggsave(karyot_time, file ="karyot_time.pdf", w = 5, h = 4)

### Generate list of individuals to sample
CM_predictions_parsed %>%
  filter(geno %in% c("Std"),
         Group == "Alys") %>%
  mutate(samp = rownames(.)  )%>% 
  as.data.frame -> Std_samps_Alys

CM_predictions_parsed %>%
  filter(geno %in% c("Std"),
         Group == "Pris") %>%
  mutate(samp = rownames(.)  )%>% 
  as.data.frame -> Std_samps_Pris

CM_predictions_parsed %>%
  filter(geno %in% c("Inv")) %>%
  mutate(samp = rownames(.)  )%>% 
  as.data.frame -> Inv_samps

#M_predictions_parsed %>%
#  filter(geno == "Het") %>%
#  mutate(samp = rownames(.)  )%>% 
#  as.data.frame -> Het_samps

write.table(Std_samps_Alys[,c("samp")], file = "Std_samps_Alys_OnlyNames.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(Std_samps_Pris[,c("samp")], file = "Std_samps_Pris_OnlyNames.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


write.table(Inv_samps[,c("samp")], file = "Inv_samps_OnlyNames.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

#write.table(Het_samps[,c("samp")], file = "Het_samps_OnlyNames.txt", 
#            quote = FALSE, sep = "\t",
#            row.names = FALSE,
#            col.names = FALSE)

#### Bring in the DGRP
dgrp_inv = fread("/project/berglandlab/DGRP_freeze2_vcf/inversion.status.txt")
dgrp_inv$`DGRP Line` = gsub("DGRP", "line", dgrp_inv$`DGRP Line`)

dgrp_inv %>%
filter(In_2L_t == "ST") %>%
  as.data.frame -> st_dgrp_samps

dgrp_inv %>%
  filter(In_2L_t == "INV") %>%
  as.data.frame -> inv_dgrp_samps

write.table(st_dgrp_samps[,c("DGRP Line")], file = "STD_DGRP_OnlyNames.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(inv_dgrp_samps[,c("DGRP Line")], file = "INV_DGRP_OnlyNames.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

##bring in the markers
#system("cp /scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt ./")

GLM_inv <- fread("./GLM_LD_outliers_annotation_priorized.txt")

data.frame(SNP_id = final_in2Lt_markers) %>%
  separate(SNP_id, into = c("chr","pos", "type")) ->
  inv_markers

rbind(cbind(inv_markers[,1:2], type = "inv"),
      cbind(GLM_inv[,c("chr","pos")], type = "glm")
      ) -> markers_of_interest


write.table(markers_of_interest, file = "Retain_loci_metadat.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(markers_of_interest, file = "Retain_loci_posOnly.txt", 
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


