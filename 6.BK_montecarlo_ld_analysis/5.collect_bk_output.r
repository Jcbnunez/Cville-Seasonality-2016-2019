### Collect the bk outfiles
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)

snp_bk_files <- system( "ls ./montecarlo_pvals_out" , intern = TRUE)

#read in files

load_in_list = list()
for(i in 1:length(snp_bk_files)){
  
  load(paste("./montecarlo_pvals_out/",
             snp_bk_files[i],
             sep = ""))
  
  out_p_df %<>%
    mutate(focal = snp_bk_files[i])
  
  load_in_list[[i]] = out_p_df
}

load_in_df = do.call(rbind, load_in_list)

load_in_df %<>%
  separate(focal,
           remove = F, 
           into = c("chr_foc", "pos_foc", "type_foc" , "test","object")) %>% 
  separate(affected_snps,
           remove = F, 
           into = c("chr_aff", "pos_aff", "type_aff" )) %>% 
  mutate(delta_pos = abs(as.numeric(pos_aff) -
                         as.numeric(pos_foc)))

load_in_df %>%
  ggplot(
    aes(
      x=delta_pos,
      y=-log10(1-p_val)
    )) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05),
             color = "red") +
  ggtitle("BK relative to inversion SNPs") ->
  test_bp

ggsave(test_bp, file = "test_bp.pdf")

####
####
####

dat = load_in_df[which(load_in_df$focal_snp == "2L_10066502_SNP"),] 

ggplot() +
  geom_violin(
    data = dat[which(dat$type == "montecarlo"),],
    color = "steelblue",
    aes(
      x=as.factor(affected_snp),
      y=pseudo_r2,
      color = type)) +
  geom_point(
    data = dat[which(dat$type == "real"),],
    fill = "red",
    shape = 5,
    aes(
      x=as.factor(affected_snp),
      y=pseudo_r2,
      color = type)) +
  coord_flip() ->
  #ggtitle(paste("Focal to", unique(tmp1$focal_snp, sep = " ") ),
  #        subtitle = paste(unique(master_obj$V1))) ->
  test_f

ggsave(test_f, file = "test_f.pdf")


####
####
####