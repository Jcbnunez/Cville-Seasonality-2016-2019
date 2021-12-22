#### Plot
#### 

library(data.table)
library(tidyverse)
library(reshape2)
library(magrittr)
library(gmodels)

load("/scratch/yey2sn/Overwintering_ms/4.GLM_snps/2L.final_output.Rdata")
chr2L =final_output
chr2L %<>% mutate(chr = "2L")

load("/scratch/yey2sn/Overwintering_ms/4.GLM_snps/2R.final_output.Rdata")
chr2R =final_output
chr2R %<>% mutate(chr = "2R")

load("/scratch/yey2sn/Overwintering_ms/4.GLM_snps/3L.final_output.Rdata")
chr3L =final_output
chr3L %<>% mutate(chr = "3L")

load("/scratch/yey2sn/Overwintering_ms/4.GLM_snps/3R.final_output.Rdata")
chr3R =final_output
chr3R %<>% mutate(chr = "3R")

rbind(chr2L,
      chr2R,
      chr3L,
      chr3R) -> final_output

final_output %>% 
  separate(Conseq, into = c("Conseq", "inv_st"), sep = "_") %>% 
  mutate(Rank = RNPV_rank,
         sets = case_when(Conseq %in% c("downstream","intergenic","Regulatory","upstream") ~ "outside_ORF",
                          !Conseq %in% c("downstream","intergenic","Regulatory","upstream") ~ "inside_ORF"),
         inv_st = case_when(inv_st == "not" ~ chr,
                            inv_st != "not" ~ inv_st)) %>% 
  group_by(Rank , Conseq, sets, inv_st, chr) %>%
  summarize(OR_m = ci(OR)[1],
            OR_l = ci(OR)[2],
            OR_h = ci(OR)[3],
  ) %>% 
  ggplot(aes(
    x=as.numeric(Rank),
    #y=log2(OR),
    y=log2(OR_m),
    ymin=log2(OR_l),
    ymax=log2(OR_h),
    color = Conseq,
    fill = Conseq
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_smooth(alpha = 0.3) +
  #geom_ribbon(alpha = 0.3) +
  geom_line() +
  geom_errorbar(width = 0.001) +
  geom_point(size=1) +
  #geom_boxplot() +
  theme_bw() +
  xlim(0, 0.11) +
  #scale_x_continuous(trans='log10') +
  facet_grid(inv_st~sets) ->
  plot_OR

ggsave(plot_OR,
       file = "plot_OR.pdf",
       width = 7,
       height = 5)


###


final_output %>% 
  separate(Conseq, into = c("Conseq", "inv_st"), sep = "_") %>% 
  mutate(Rank = RNPV_rank,
         sets = case_when(Conseq %in% c("downstream","intergenic","Regulatory","upstream") ~ "outside_ORF",
                          !Conseq %in% c("downstream","intergenic","Regulatory","upstream") ~ "inside_ORF",
         )) %>%
  ggplot(aes(
    x=as.numeric(Rank),
    y=-log10(p.val),
    #y=log2(OR_m),
    color = Conseq,
    fill = Conseq,
    linetype = chr
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #geom_smooth(alpha = 0.3) +
  #geom_ribbon(alpha = 0.3) +
  geom_line() +
  #geom_errorbar(width = 0.001) +
  #geom_point(size=1) +
  #geom_boxplot() +
  theme_bw() +
  facet_grid(inv_st~sets) ->
  plot_OR

ggsave(plot_OR,
       file = "plot_OR.box.pdf",
       width = 6,
       height = 4)

