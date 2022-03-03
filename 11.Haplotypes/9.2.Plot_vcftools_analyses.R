### Collect VCFTools analyses
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)

inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"
inv.dt <- fread(inversion_map)
setnames(inv.dt, "chrom", "chr.x")

####### Load all

load("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/pi_D_datforplot.Rdata")

all_pi_d_parsed %>%
  ggplot(
    aes(
      x=mid_point,
      y=value,
      #color =karyo,
      linetype = pop)
  ) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_line(alpha = 0.5) +
  ggtitle("All Samples") +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  pi_d_plot_all

ggsave(pi_d_plot_all, file = "pi_d_plot_all.pdf")

sub_pi_d_parsed %>% 
  group_by(pop, type, resolution) %>%
  summarize(N = n())

sub_pi_d_parsed %>%
  filter(resolution == "W_100000") %>%
  ggplot(
    aes(
      x=BIN_START,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=start), linetype="dashed") +
  geom_vline(data=inv.dt[which(inv.dt$invName == "2Lt")], aes(xintercept=stop), linetype="dashed") +
  geom_line(alpha = 0.5) +
  ggtitle("Samples split by Karyotype") +
  theme_bw() +
  facet_wrap(~variable, ncol = 1, scales = "free") ->
  pi_d_plot_sub

ggsave(pi_d_plot_sub, file = "pi_d_plot_sub.pdf")

###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
###### zoom Window
#####  ######
load("/scratch/yey2sn/Overwintering_ms/11.Haplotypes/pi_D_datforplot.Rdata")


sub_pi_d_parsed %>%
  filter(resolution == "W_10000") %>% 
  mutate(zoom_win = case_when(BIN_START > 4.7e6 & BIN_START < 5.3e6 ~ "1.win5",
                              BIN_START > 5.8e6 & BIN_START < 7e6 ~ "2.win6",
                              BIN_START > 9.0e6 & BIN_START < 10.6e6 ~ "3.win10",
                              )) %>% 
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) -> D_data_zoom

D_data_zoom %>%
  ggplot(
      aes(
      x=BIN_START/1e6,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_line(alpha = 0.5) +
  ggtitle("Zoom-in to GLM peaks") +
  theme_bw() +
  facet_wrap(variable~zoom_win, ncol = 3, scales = "free") ->
  pi_d_plot_zoom

ggsave(pi_d_plot_zoom, file = "pi_d_plot_zoom.pdf", h= 6, w = 9)


###
### add genes regions
#bring in the snp data:
glm_ld_outliers <- "/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt"
outlier_list <- fread(glm_ld_outliers)

outlier_list %>%
  filter(Transcript_biotype == "protein_coding") %>%
  group_by(Gene_Name) %>%
  summarise(min_pos = min(pos),
            max_pos = max(pos))  %>%
  mutate(zoom_win = 
           case_when(min_pos > 4.7e6 & max_pos < 5.3e6 ~ "1.win5",
                     min_pos > 5.8e6 & max_pos < 7e6 ~ "2.win6",
                     min_pos > 9.0e6 & max_pos < 10.6e6 ~ "3.win10")) -> gene_list

gene_list  %>% 
  .[-grep("^CG" ,.$Gene_Name),] %>%
  filter(zoom_win %in% c("1.win5", "2.win6","3.win10" )) %>%
  .[order(.$min_pos),] %>%
  mutate(gene_id = 1:dim(.)[1],
         analysis = "2.GeneId") ->
  gene_list_select

D_data_zoom %<>%
  mutate(analysis = "1.TajD")

data.frame(xint = c(6.07, 6.69, 9.64, 10.0),
           zoom_win = c("2.win6", "2.win6", "3.win10","3.win10")) ->
  vlines

ggplot() +
  geom_vline(data=vlines, aes(xintercept=xint), linetype="dashed") +
  geom_line(data= filter(D_data_zoom, variable == "TajimaD", pop == "cm"),
    aes(
      x=BIN_START/1e6,
      y=value,
      color =type,
      linetype = pop)
  ) + 
  geom_segment(data=gene_list_select,
               aes(
                x=min_pos/1e6,
                xend=max_pos/1e6,
                y=gene_id/1e6,
                yend=gene_id/1e6), alpha = 0.5
               )  +
  geom_text(data=gene_list_select,
               aes(
                 x=(max_pos)/1e6,
                 y=(gene_id+0.5)/1e6,
                 label=Gene_Name
                 ),size = 1
  )  +
  facet_wrap(analysis~zoom_win, ncol = 3, scales = "free") ->
  gene_d_plot_zoom

ggsave(gene_d_plot_zoom, file = "gene_d_plot_zoom.pdf", h= 6, w = 9)

  

