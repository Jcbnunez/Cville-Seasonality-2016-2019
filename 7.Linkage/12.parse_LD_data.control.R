### Code to process ld data
### 
library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)
library(car)
library(DescTools)
#parallel computing in R
library(foreach)
library(doMC)
registerDoMC(10) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem
###
##########
##########
####### ------> can start in line 55
####### 
## load the parsed object
load("./MATCHEDCONTROLs.merged.ld.Rdata")

### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B))


#save(matched_controls, GLM_LD_metadata_object, file = "macthed_controls_interm.file.Rdata")
load("./macthed_controls_interm.file.Rdata")
#matched_controls
#GLM_LD_metadata_object
ld_df %>%
  mutate(type_A = case_when(BP_A %in% GLM_LD_metadata_object$pos ~ "control.A")) %>%
  mutate(type_B = case_when(BP_B %in% GLM_LD_metadata_object$pos ~ "control.B")) ->
    ld_df_typ

ld_df_typ$type_A[is.na(ld_df_typ$type_A)] = "control"
ld_df_typ$type_B[is.na(ld_df_typ$type_B)] = "control"

ld_df_typ %<>% mutate(match =  paste(type_A, type_B, sep = "_") )
ld_df_typ_control = ld_df_typ
######

#### LD plot
load("./merged.ld.Rdata")
### Add the bp distance between SNPs
ld_df %<>% 
  mutate(BP_diff = abs(BP_A-BP_B))

ld_df %>%
  filter(BP_A %in% GLM_LD_metadata_object$pos & BP_B %in% GLM_LD_metadata_object$pos ) %>% 
  mutate(match = "GLM_GLM") ->
  ld_df_glm

rbind(ld_df_typ_control, ld_df_glm, fill = T) ->
  ld_dat_merged
############################ start here
########################################################
############################
############################
############################
############################
###
#save(ld_dat_merged, file ="ld_dat_glm_control_merged.Rdata")

load("./ld_dat_glm_control_merged.Rdata")
### ---> loads ld_dat_merged

ld_dat_merged %>%
mutate(dist_bin = RoundTo(BP_diff, 100000, "floor"),
       BPA_binned = RoundTo(BP_A, 100000, "floor"),
       ) -> 
  ld_dat_merged_binned


###
ld_dat_merged_binned %>%
  group_by(BPA_binned, match) %>%
  summarize(Mean_pos = mean(R2),
            SD = sd(R2)) ->
  ld_dat_merged_binned_mean_per_pos

####
ld_dat_merged_binned_mean_per_pos %>%
  ggplot(
    aes(x=BPA_binned,
        y=Mean_pos,
        ymin=Mean_pos+SD,
        ymax=Mean_pos-SD,
        color=match,
        fill=match)
  ) + 
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  geom_point() +
  ylab("Mean R2 of each 50k window") +
  xlab("Window Start") ->
  ld_plots_means

ggsave(ld_plots_means, file = "ld_plots_means.pdf")


#####
ld_dat_merged_binned %>% 
  group_by( match,  dist_bin) %>%
  summarize(N=n(),
            Mean = mean(R2),
            SD = sd(R2)) ->
  ld_summarized

######
ld_summarized %>%
  ggplot(
    aes(
      x=dist_bin,
      y=Mean,
      ymin=Mean+SD,
      ymax=Mean-SD,
      fill=match,
      color=match)
  ) +
  geom_ribbon(alpha = 0.5) +
  geom_point() ->
  ld_plot_smooth

ggsave(ld_plot_smooth, file = "ld_plot_smooth.pdf")





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



##################### Analyzes start here
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


### find long distance pairs
### 
### 
### 
### General parameters
R2_cutoff = 0.6
low_distance = 1e6
#high_diatnce = 9e6

### Make LD filtered object
ld_df %>%
  filter(Comp_simplified == "vs_out" & BP_diff > low_distance | Comp_simplified == "vs_inv",
         R2 > R2_cutoff) ->
  LD_test_out

#save SNPs pairs of interest
save(LD_test_out, file = "ld_df_intermediate_step.Rdata")

### see individual SNPs
LD_test_out %>%
  .$BP_A -> bp_a

LD_test_out %>%
  .$BP_B -> bp_b

paste( "there are", length(unique(c(bp_a, bp_b))), " SNPs that pass filter" )

## create data frame with long distance partner outliers
data.frame(chr = "2L",
           pos = unique(c(bp_a, bp_b)) ) ->
  long_distance_partners

inv_markers_id %>%
  separate(V1, into = c("chr", "pos", "feature"), sep = "_") %>%
  mutate(pos = as.numeric(.$pos),
         inversion_marker = "TRUE") %>%
  dplyr::select(chr, pos, inversion_marker) %>%
  right_join(long_distance_partners) ->
  long_distance_partners_inv


## Save table for GOWINDA
write.table(long_distance_partners_inv,
            file = paste(R2_cutoff,"GLM_LD_outliers.txt", sep = "_"), 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = TRUE, 
            qmethod = c("escape", "double"),
            fileEncoding = "")



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







