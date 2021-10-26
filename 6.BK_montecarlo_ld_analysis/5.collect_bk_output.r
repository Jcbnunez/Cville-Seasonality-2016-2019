### Collect the bk outfiles
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)

snp_bk_files <- system( "ls ./bk_out_sims" , intern = TRUE)

#read in files

load_in_list = list()
for(i in 1:length(snp_bk_files)){
  
  load(paste("./bk_out_sims/",
             snp_bk_files[i],
             sep = ""))
  
  load_in_list[[i]] = results_df
}

load_in_df = do.call(rbind, load_in_list)


# Estimate p-values:
dat_in = load_in_df

dat_in %>%
  filter(test_type == "inv>inv",
         type == "real") -> real_inv

dat_in %>%
  filter(test_type == "inv>inv",
         type == "montecarlo") -> mont_inv

ggplot() +
  geom_boxplot(
    data = mont_inv,
    aes(x= affected_snp,
        y= pseudo_r2),
    outlier.shape = NA
  ) +
  geom_point(
    data = real_inv,
    aes(x=affected_snp,
        y= pseudo_r2),
    size = 3,
    color = "red",
    shape = 18
    
  ) +
  coord_flip() +
  ggtitle("Do Inversion predict themselves?") +
  ylim(0,0.4) +
  facet_wrap(~focal_snp) ->
  test

ggsave(test,
       file = "inversions_bby_inversion.dist.png"
       )

ggplot(mont_inv, aes(x=x, y=y) ) +
  geom_hex() +
  theme_bw()


######### selection gml
# add link to data
obj <- "/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata"
#load the object
load(obj)

## extarct outlier SNPs
glm.out %>%
  filter(mod=="aveTemp",
         rnp.clean<0.001,
         chr=="2L") ->
  glm.out_p1

glm.out_p1 %<>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))



dat_in %>%
  filter(test_type == "inv>tmp",
         type == "real") -> real_tmp

dat_in %>%
  filter(test_type == "inv>tmp",
         type == "montecarlo") -> mont_tmp

real_tmp$focal_snp %>% unique -> lists1
real_tmp$affected_snp %>% unique -> lists2

ggplot() +
  geom_boxplot(
    data = mont_tmp[which(mont_tmp$affected_snp %in% glm.out_p1$SNP_id),],
    aes(x= affected_snp,
        y= pseudo_r2),
    outlier.shape = NA
  ) +
  geom_point(
    data = real_tmp[which(real_tmp$affected_snp %in% glm.out_p1$SNP_id),],
    aes(x=affected_snp,
        y= pseudo_r2),
    size = 3,
    color = "red",
    shape = 18
  ) +
  ggtitle("Do Inversion predict glm?") +
  coord_flip() +
  facet_wrap(~focal_snp) ->
  test2

ggsave(test2,
       file = "inversions_by_glm.dist.png")

##### Calculate P-values


########
output_list = list()
for(i in 1:length(focal_snps)){
  
  print(i/length(focal_snps)*100)  
  dat_in %>%
    filter(focal_snp == focal_snps[i]) -> tmp
  
  tmp %>% 
    filter(type == "real") %>%
    group_by(affected_snp, test_type) %>%
    summarize(Real_R = mean(pseudo_r2)) %>% 
    mutate(focal_snp = focal_snps[i])->
    tmp_real
  
  tmp_real$prediction_robustness = NA
  #tmp_real$t.test.mu = NA
  
  for(j in 1:dim(tmp_real)[1]){
    
    tmp_real$prediction_robustness[j] = mean(tmp[which(tmp$affected_snp == tmp_real$affected_snp[j] & 
                                                         tmp$type == "montecarlo" ),]$pseudo_r2 >= 
                                               tmp_real$Real_R[j])

  } ## close j
  
  #tmp_real %<>%
  #  mutate(reversed_t = 1 - as.numeric(t.test.mu))
  
  output_list[[i]] = tmp_real
} ## close i

######
bk_output = do.call(rbind, output_list )

save(load_in_df, bk_output,
     file = "bk.test.output.Rdata")


