
#### Bootstrap and resampling analysis


#load packages
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)
library(gtools)
library(patchwork)
library(ggalt)
library(reshape2)
library(ggrepel)
library(sp)
library(permute)
library(RColorBrewer)
library(devtools)
library(lubridate)
#install_github('tavareshugo/windowscanr')
library(windowscanr)

load("./var_fst_calc.Rdata")
variance_of_FST %>% head
names(variance_of_FST)[1] = "Month"

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

sample_n =  SNP_sampler[k]   
  
loop1_list = list()
for(i in 1:length(objects)){
  
  load(objects[i])
  o %>% colnames() %>% .[1] -> lead_snp
  
  lead_snp %>% data.frame(headsnp = .) %>% 
    separate(headsnp , into = c("CHR","POS")) -> lead_snp_guide
  
  print(lead_snp_guide$CHR)
  print(sample_n)
  
  filtered_samps_for_analysis %>%
    filter(city == "Charlottesville",
           MeanEC > 30) %>% 
    .$sampleId -> select_samples
  
  o %>%
    as.data.frame() %>% 
    filter(rownames(.) %in%  select_samples) %>% 
    .[,sample(dim(o)[2],sample_n )] %>% 
    PCA(scale.unit = F, graph = F, ncp = 20) ->
    PCA_object
  
  PCA_object$ind$coord %>%
    as.data.frame() %>%
    mutate(sampleId = rownames(.),
           chr = lead_snp_guide$CHR,
           snp_n = sample_n) %>%  
    left_join(., filtered_samps_for_analysis ) %>%
    left_join(variance_of_FST) -> PCA_table
  
  ### Battery of tests
rbind(
  data.frame(
    pc=1,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="year",
    cor=cor.test(PCA_table$Dim.1, PCA_table$year)$estimate^2
    ),
  data.frame(
    pc=2,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="year",
    cor=cor.test(PCA_table$Dim.2, PCA_table$year)$estimate^2
  ),
  data.frame(
    pc=3,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="year",
    cor=cor.test(PCA_table$Dim.3, PCA_table$year)$estimate^2
  )
) -> time_correlation

##inversion
rbind(
  data.frame(
    pc=1,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="inv",
    cor=cor.test(PCA_table$Dim.1, PCA_table$`In(2L)t`)$estimate^2
  ),
  data.frame(
    pc=2,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="inv",
    cor=cor.test(PCA_table$Dim.2, PCA_table$`In(2L)t`)$estimate^2
  ),
  data.frame(
    pc=3,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="inv",
    cor=cor.test(PCA_table$Dim.3, PCA_table$`In(2L)t`)$estimate^2
  )
) -> inv_correlation

#Nc
rbind(
  data.frame(
    pc=1,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="Nc",
    cor=cor.test(PCA_table$Dim.1, PCA_table$MeanEC)$estimate^2
  ),
  data.frame(
    pc=2,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="Nc",
    cor=cor.test(PCA_table$Dim.2, PCA_table$MeanEC)$estimate^2
  ),
  data.frame(
    pc=3,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="Nc",
    cor=cor.test(PCA_table$Dim.3, PCA_table$MeanEC)$estimate^2
  )
) -> MeanEC_correlation

#lm - month
rbind(
  data.frame(
    pc=1,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="lm_month_2",
    cor=summary(lm(abs(Dim.1) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6:12),] ) )$adj.r.squared
  ),
  data.frame(
    pc=2,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="lm_month_2",
    cor=summary(lm(abs(Dim.2) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6:12),] ) )$adj.r.squared
  ),
  data.frame(
    pc=3,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="lm_month_2",
    cor=summary(lm(abs(Dim.3) ~ Month + I(Month^2), data = PCA_table[which(PCA_table$chr == lead_snp_guide$CHR & PCA_table$Month %in% 6:12),] ) )$adj.r.squared
  )
) -> month_correlation

##Variance_object 
PCA_table %>%
  select(FST_sd, 
         Month, 
         Dim.1,
         Dim.2,
         Dim.3) %>%
  melt(id = c("FST_sd","Month")) %>% 
  .[complete.cases(.),] %>% 
  group_by(FST_sd, Month, variable) %>%
  summarise(PC_var = sd(value)) ->
  variance_object

rbind(
  data.frame(
    pc=1,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="VarFst_Var_PC",
    cor=cor.test(variance_object$FST_sd[which(variance_object$variable == "Dim.1")] , 
                 variance_object$PC_var[which(variance_object$variable == "Dim.1")])$estimate
  ),
  data.frame(
    pc=2,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="VarFst_Var_PC",
    cor=cor.test(variance_object$FST_sd[which(variance_object$variable == "Dim.2")] , 
                 variance_object$PC_var[which(variance_object$variable == "Dim.2")])$estimate
  ),
  data.frame(
    pc=3,
    chr= lead_snp_guide$CHR,
    snp_n = sample_n,
    test="VarFst_Var_PC",
    cor=cor.test(variance_object$FST_sd[which(variance_object$variable == "Dim.3")] , 
                 variance_object$PC_var[which(variance_object$variable == "Dim.3")])$estimate
  )
) -> Var_correlation

### joint all objects

rbind(
time_correlation,
inv_correlation,
MeanEC_correlation,
month_correlation,
Var_correlation
) %>%
  as.data.frame() %>% 
  mutate(Type = "Real") ->
  tmp_final

loop1_list[[i]] = tmp_final
} # i

loop2_list[[k]] = do.call(rbind, tmp_final)

} # k




######
######


pca_table_df %>%
  filter(year >= 2016) %>%
  group_by(month_col, chr) %>%
  summarise(sum_dim1 = sd(Dim.1),
            sum_dim2 = sd(Dim.2),
            sum_dim3 = sd(Dim.3)) %>% 
  left_join(variance_of_FST) ->
  var_per_month

var_per_month %>%
  filter(month_col %in% 7:11) %>% 
  melt(id = c("month_col", "chr", "FST_sd")) %>% 
  ggplot(aes(
    x=FST_sd,
    y=value,
    shape = chr,
    color = month_col
  )) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", 
              color = "black",
              se = F) +
  xlab(expression(Var(italic(F)[ST]))) +
  ylab("Var(PC projection)") +
  facet_grid(chr~variable) ->
  var_month

ggsave(var_month,
       file = "var_month.pdf")

var_per_month %>%
  filter(year > 2016,
         month_col %in% 7:11) ->
  data_in

lm(sum_dim2 ~ month_col + I(month_col^2), data = data_in[which(data_in$chr == "3R"),] ) ->
  lm_sq_month
summary(lm_sq_month)


