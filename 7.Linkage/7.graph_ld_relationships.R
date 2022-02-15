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

## import the ld object
load("./ld_df_intermediate_step.Rdata")
# -- > LD_test_out

## import functional data
func_prioritized <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

#### prepare join
func_prioritized %>%
  dplyr::select(CHR_A = chr, BP_A = pos, Annot_A = Annotation, Gene_A = Gene_Name) -> A_tab

func_prioritized %>%
  dplyr::select(CHR_B = chr, BP_B = pos, Annot_B = Annotation, Gene_B = Gene_Name) -> B_tab

### Merge
left_join(LD_test_out, A_tab) %>% 
  left_join(B_tab) %>%
  mutate(Joint_annot = paste(Annot_A, Annot_B, sep = "|")) -> 
  LD_test_out_annot

LD_test_out_annot$Joint_annot = gsub("_variant", "", LD_test_out_annot$Joint_annot )


### long range targets
###       pos >= 2225744, 
#         pos <= 13154180

LD_test_out_annot %<>%
  mutate(rank = order(LD_test_out$BP_A))

LD_test_out_annot %>%
  filter(
    Annot_A == "missense_variant" & Annot_B == "missense_variant",
    BP_A < BP_B,
    R2 > 0.6,
    BP_A > 3e6,
    BP_A < 12e6,
    BP_B > 3e6,
    BP_B < 12e6,) -> missense_pairs

missense_pairs %>%
  filter(Gene_A == "mtDNA-helicase" | Gene_B == "mtDNA-helicase" ) %>%
  mutate(gene_pair = case_when(Gene_A == "mtDNA-helicase" ~ Gene_B,
                               Gene_B == "mtDNA-helicase" ~ Gene_A)  ) %>%
  ggplot(aes(
    x = BP_A,
    xend = BP_B,
    y= 1:dim(.)[1],
    yend= 1:dim(.)[1],
    label = gene_pair,
    color = R2
  )) +
  ggtitle("SNPs in strong linkage (LD > 0.6) to the Mitochondrial helicase gene",
          subtitle = "Dashed line shows mtDNA-helicase") +
  geom_segment(size = 2) +
  geom_text(size = 3, nudge_x = -700000, color = "black") +
  xlim(0,15e6) + 
  geom_vline(xintercept = (9891255+9893446)/2, linetype = "dashed") +
  scale_color_gradient2(low = "purple", high = "red", mid = "gold", midpoint =  0.70) +
  #facet_wrap(~Gene_A, scales = "free_y") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  LD_ranges_filter_mito

ggsave(LD_ranges_filter_mito, 
       file = "LD_ranges_filter_mito.png",
       width = 7,
       height = 7)

missense_pairs %>%
  ggplot(aes(
    x = BP_A,
    xend = BP_B,
    y= 1:dim(.)[1],
    yend= 1:dim(.)[1],
    label = Gene_B,
    color = R2
  )) +
  ggtitle("SNPs in strong linkage (LD > 0.6) that are non-synonymous") +
  geom_segment(size = 1) +
  geom_text(size = 1, nudge_x = -700000, color = "black") +
  xlim(0,15e6) + 
  scale_color_gradient2(low = "purple", high = "red", mid = "gold", midpoint =  0.70) +
  facet_wrap(~Gene_A, scales = "free") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  LD_ranges_filter

ggsave(LD_ranges_filter, 
       file = "LD_ranges_filter.png",
       width = 15,
       height = 15)