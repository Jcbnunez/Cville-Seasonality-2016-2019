
#### Bootstrap and resampling analysis
rm(list = ls())

#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
#library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
#install_github('tavareshugo/windowscanr')
#library(windowscanr)

args = commandArgs(trailingOnly=TRUE)

array_num=args[1]

root="/scratch/yey2sn/Overwintering_ms/3.Cville_PCA/real_resample_out"

objects <- c(
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2R.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3L.ECfiltered.Rdata",
  "/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_3R.ECfiltered.Rdata"
)

SNP_sampler <- c(seq(from = 100, to = 1000, by = 100),
                 seq(from = 1000, to = 20000, by = 1000))


# run loop
loop2_list = list()
for(k in 1:length(SNP_sampler)){
  
  sample_nsnps =  SNP_sampler[k]   
  
  loop1_list = list()
  for(i in 1:length(objects)){
    
    load(objects[i])
    o %>% colnames() %>% .[1] -> lead_snp
    
    lead_snp %>% data.frame(headsnp = .) %>% 
      separate(headsnp , into = c("CHR","POS")) -> lead_snp_guide
    
    print(lead_snp_guide$CHR)
    print(sample_nsnps)
    
    filtered_samps_for_analysis %>%
      filter(city == "Charlottesville",
             MeanEC > 30) %>% 
      .$sampleId -> select_samples
    
    o %>%
      as.data.frame() %>% 
      filter(rownames(.) %in%  select_samples) %>% 
      .[,sample(dim(o)[2],sample_nsnps )] %>% 
      PCA(scale.unit = F, graph = F, ncp = 20) ->
      PCA_object
    
    PCA_object$ind$coord %>%
      as.data.frame() %>%
      mutate(sampleId = rownames(.),
             chr = lead_snp_guide$CHR,
             snp_n = sample_nsnps) %>%  
      left_join(., filtered_samps_for_analysis ) -> PCA_table
    
    ### Battery of tests
    rbind(
      data.frame(
        pc=1,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="year",
        cor=cor.test(PCA_table$Dim.1, PCA_table$year)$estimate,
        cor_low=cor.test(PCA_table$Dim.1, PCA_table$year)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.1, PCA_table$year)$conf.int[2]
      ),
      data.frame(
        pc=2,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="year",
        cor=cor.test(PCA_table$Dim.2, PCA_table$year)$estimate,
        cor_low=cor.test(PCA_table$Dim.2, PCA_table$year)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.2, PCA_table$year)$conf.int[2]
      ),
      data.frame(
        pc=3,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="year",
        cor=cor.test(PCA_table$Dim.3, PCA_table$year)$estimate,
        cor_low=cor.test(PCA_table$Dim.3, PCA_table$year)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.3, PCA_table$year)$conf.int[2]
      )
    ) -> time_correlation
    
    ##inversion
    
    rbind(
      data.frame(
        pc=1,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="In(2L)t",
        cor=cor.test(PCA_table$Dim.1, PCA_table$`In(2L)t`)$estimate,
        cor_low=cor.test(PCA_table$Dim.1, PCA_table$`In(2L)t`)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.1, PCA_table$`In(2L)t`)$conf.int[2]
      ),
      data.frame(
        pc=2,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="In(2L)t",
        cor=cor.test(PCA_table$Dim.2, PCA_table$`In(2L)t`)$estimate,
        cor_low=cor.test(PCA_table$Dim.2, PCA_table$`In(2L)t`)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.2, PCA_table$`In(2L)t`)$conf.int[2]
      ),
      data.frame(
        pc=3,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="In(2L)t",
        cor=cor.test(PCA_table$Dim.3, PCA_table$`In(2L)t`)$estimate,
        cor_low=cor.test(PCA_table$Dim.3, PCA_table$`In(2L)t`)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.3, PCA_table$`In(2L)t`)$conf.int[2]
      )
    ) -> inv_correlation
    
    
    #Nc
    rbind(
      data.frame(
        pc=1,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="EFC",
        cor=cor.test(PCA_table$Dim.1, PCA_table$MeanEC)$estimate,
        cor_low=cor.test(PCA_table$Dim.1, PCA_table$MeanEC)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.1, PCA_table$MeanEC)$conf.int[2]
      ),
      data.frame(
        pc=2,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="EFC",
        cor=cor.test(PCA_table$Dim.2, PCA_table$MeanEC)$estimate,
        cor_low=cor.test(PCA_table$Dim.2, PCA_table$MeanEC)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.2, PCA_table$MeanEC)$conf.int[2]
      ),
      data.frame(
        pc=3,
        chr= lead_snp_guide$CHR,
        snp_n = sample_nsnps,
        test="EFC",
        cor=cor.test(PCA_table$Dim.3, PCA_table$MeanEC)$estimate,
        cor_low=cor.test(PCA_table$Dim.3, PCA_table$MeanEC)$conf.int[1],
        cor_high=cor.test(PCA_table$Dim.3, PCA_table$MeanEC)$conf.int[2]
      )
    ) -> MeanEC_correlation
    
    ###lm - month
    ##rbind(
    ##  data.frame(
    ##    pc=1,
    ##    chr= lead_snp_guide$CHR,
    ##    snp_n = sample_nsnps,
    ##    test="lm_month_2",
    ##    cor=summary(lm(abs(Dim.1) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6##:12),] ) )$r.squared,
    ##    cor_low=summary(lm(abs(Dim.1) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared,
    ##    cor_high=summary(lm(abs(Dim.1) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared
    ##  ),
    ##  data.frame(
    ##    pc=2,
    ##    chr= lead_snp_guide$CHR,
    ##    snp_n = sample_nsnps,
    ##    test="lm_month_2",
    ##    cor=summary(lm(abs(Dim.2) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6##:12),] ) )$r.squared,
    ##    cor_low=summary(lm(abs(Dim.2) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared,
    ##    cor_high=summary(lm(abs(Dim.2) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared
    ##  ),
    ##  data.frame(
    ##    pc=3,
    ##    chr= lead_snp_guide$CHR,
    ##    snp_n = sample_nsnps,
    ##    test="lm_month_2",
    ##    cor=summary(lm(abs(Dim.3) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6##:12),] ) )$r.squared,
    ##    cor_low=summary(lm(abs(Dim.3) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared,
    ##    cor_high=summary(lm(abs(Dim.3) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month ##%in% 6:12),] ) )$r.squared
    ##  )
    ##) -> month_correlation
    
    ### joint all objects
    
    rbind(
      time_correlation,
      inv_correlation,
      MeanEC_correlation,
      ##month_correlation
    ) %>%
      as.data.frame() %>% 
      mutate(Type = "Real",
             run_id = array_num) ->
      tmp_final
    
    loop1_list[[i]] = tmp_final
  } # i
  
  loop2_list[[k]] = do.call(rbind, loop1_list)
  
} # k

real_df = do.call(rbind, loop2_list)

write.table(real_df, 
            file = paste(root,"/real.resample.",array_num,".txt", sep =""),
            append = FALSE,
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = TRUE)

