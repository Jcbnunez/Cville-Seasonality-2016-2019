### Code to process ld data
### 

library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)
#parallel computing in R
library(foreach)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem
###
##########

## load the parsed object
load("./merged.ld.Rdata")

### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B))

## bring in the inversion data -- for labeling SNPS.
inv_markers_id <- fread("/scratch/yey2sn/Overwintering_ms/Inversion_markers/in2lt_ld_47snps_informative_markers.txt", head = F)

## Add inv/outlier label
ld_df %<>%
  mutate(A_status = case_when(SNP_A %in% inv_markers_id$V1 ~ "inv",
                              !(SNP_A %in% inv_markers_id$V1) ~ "out"),
         B_status = case_when(SNP_B %in% inv_markers_id$V1 ~ "inv",
                              !(SNP_B %in% inv_markers_id$V1) ~ "out"))  %>% 
  mutate(Comp_stat = paste(A_status, B_status, sep = "_")) %>%
  mutate(Comp_simplified = case_when(A_status == "inv" | B_status == "inv" ~ "vs_inv",
                                     A_status != "inv" & B_status != "inv" ~ "vs_out"))


##sanity checks
ld_df %>%
  group_by(Comp_stat, Comp_simplified) %>%
  summarise(N=n())

### 
### find the distribution of LD
ld_df %>%
  ggplot(aes(R2)) +
  geom_histogram() ->
  ld_histogram

ggsave(ld_histogram, file = "ld_histogram.pdf")

#####
ld_df %>%
  mutate(BP_diff_scaled = round(BP_diff/1e6, 1)) %>% 
  group_by(BP_diff_scaled, Comp_simplified) %>%
  summarise(mean_ld = mean(R2),
            sd_ld = sd(R2),
            median_ld = quantile(R2, 0.5), 
            ld_lowQ = quantile(R2, 0.05), 
            ld_highQ = quantile(R2, 0.95)) ->
  summaries_of_ld

summaries_of_ld %>%
  ggplot(aes(
    x=BP_diff_scaled,
    y=median_ld,
    ymin =  ld_lowQ,
    ymax =  ld_highQ
  )) + 
  geom_ribbon(alpha = 0.4) + 
  geom_line() +
  facet_grid(~Comp_simplified)->
  ld_summaries_fig

ggsave(ld_summaries_fig, file = "ld_summaries_fig.pdf")

#### Plot LD triangle
ld_df %>%
  filter(R2 > 0.1,
         BP_A < BP_B) -> all_dat

ld_df %>%
  filter(R2 > 0.7,
         BP_A < BP_B) -> high_dat


### This is an output from code #6... may need tor un that first
func_prioritized <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

left_join(func_prioritized, 
          data.frame(pos=unique(high_dat$BP_A))) %>%
  group_by(Gene_Name) %>%
  summarise(mean_pos = mean(pos)) -> a

  hist(a$mean_pos, breaks = 30) -> b
data.frame(b$density, b$mids) %>%
  ggplot(aes(x=b.mids, y=b.density)) +
  geom_line() -> ld_pairs

ggsave(ld_pairs, file = "ld_pairs.pdf")



  ggplot() +
    geom_point(data=all_dat,
                 aes(x=BP_A,
                   y=BP_B,
                   color = R2),
               size = 1, alpha = 0.5, shape = 15) +
    geom_point(data=high_dat,
               aes(x=BP_A,
                   y=BP_B,
                   color = R2),
               size = 0.1, alpha = 0.5, shape = 15) +
    geom_vline(xintercept = 5500000, linetype = "dashed") +
    geom_vline(xintercept = 6500000, linetype = "dashed") +
    geom_vline(xintercept = 9500000, linetype = "dashed") +
    theme_classic() +
  scale_color_gradient2(low = "blue", high = "red", mid = "gold", midpoint =  0.5) ->
  ld_triag_plot

ggsave(ld_triag_plot, file = "raster_plot.png")


### find long distance pairs
### 
### 
### 
### General parameters
R2_cutoff = 0.6
low_distance = 1e6
high_diatnce = 9e6


### Make LD filtered object
ld_df %>%
  filter(Comp_simplified == "vs_out",
         BP_diff > low_distance,
         BP_diff < high_diatnce,
         R2 > R2_cutoff) ->
  LD_test_out

save(LD_test_out, file = "ld_df_intermediate_step.Rdata")




####

LD_test_out %>%
  ggplot(aes(x=BP_diff,
             y=R2)) +
  geom_point() ->
  long_distance_ld_plot

ggsave(long_distance_ld_plot, file = "long_distance_ld_plot.png")

### See distribution of delta bp
LD_test_out %>%
  ggplot(aes(BP_diff)) +
  geom_histogram() ->
  hist_bd_diff_filt

ggsave(hist_bd_diff_filt, file = "hist_bd_diff_filt.png")

### see individual SNPs
LD_test_out %>%
  .$BP_A -> bp_a

LD_test_out %>%
  .$BP_B -> bp_b

paste( "there are", length(unique(c(bp_a, bp_b))), " SNPs that pass filter" )

## create data frame with long distance partner outliers
data.frame(chr = "2L",
           pos = unique(c(bp_a, bp_b)) ) %>%
  long_distance_partners

## Save table for GOWINDA
write.table(long_distance_partners,
            file = "./GLM_LD_outliers.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")





