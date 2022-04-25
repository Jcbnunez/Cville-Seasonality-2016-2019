### Extarct SNPs within windows of interest
### 

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(ggrepel)

annots_priority <- fread("/scratch/yey2sn/Overwintering_ms/7.LD/GLM_LD_outliers_annotation_priorized.txt")

annots_priority %<>%
mutate(window_of_int = 
         case_when(pos > 5155762 & pos < 5255762 ~ "1.win5",
                   pos > 6255762 & pos < 6355762 ~ "2.win6",
                   pos > 9505762 & pos < 9605762 ~ "3.win9")) 


annots_priority %>% 
  group_by(window_of_int) %>%
  summarize(N = n()) %>% 
  as.data.frame() %>%
  .[complete.cases(.),]

annots_priority %>% 
  group_by(window_of_int, Annotation) %>%
  summarize(N = n()) %>% 
  as.data.frame() %>%
  .[complete.cases(.),]

annots_priority %>% 
  filter(window_of_int %in% c("1.win5", "2.win6", "3.win9"),
         #Annotation == "missense_variant"
  ) %>%
  group_by(Gene_Name) %>%
  slice_head() %>% arrange(pos) %>%   as.data.frame() 


annots_priority %>% 
  filter(window_of_int %in% c("1.win5", "2.win6", "3.win9"),
         Annotation == "missense_variant"
         ) %>%
  group_by(Gene_Name) %>%
  slice_head() %>% arrange(pos) %>%   as.data.frame() 

##genes of interest

#Msp300
#CG14034 -- Triacylglycerol lipase 

#ppk14 -- 

#CG3841 -- Carboxylesterase

### see raw p-values

glm.file <- "/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_0.Rdata"

load(glm.file)

glm.out %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_")) -> 
  glm.ids

glm.ids %>%
  filter(mod == "aveTemp+year_factor") %>%
  dim

adjst=50000
glm.ids %>%
  filter(chr == "2L",
         pos > 5e6 & pos < 10e6,
         mod == "aveTemp+year_factor") %>% 
  mutate(window_of_int = 
           case_when(pos > 5155762-adjst & pos < 5255762+adjst ~ "1.win5",
                     pos > 6255762-adjst & pos < 6355762+adjst ~ "2.win6",
                     pos > 9505762-adjst & pos < 9605762+adjst ~ "3.win9")) ->
  glm.ids.winInts

starts_bp = c(5155762, 6255762, 9505762)
ends_bp = c(5255762, 6355762, 9605762)

starts_snp = c(5144479, 6309567, 9568122)
ends_snp = c(5218310, 6356630, 9611648)


left_join(
  filter(annots_priority, Annotation %in% c("synonymous_variant","5_prime_UTR_variant","3_prime_UTR_variant","missense_variant")), 
  filter(glm.ids.winInts, window_of_int %in% c("1.win5", "2.win6", "3.win9") )) -> outs

### plot
  ggplot() +
  geom_point(data = filter(glm.ids.winInts, window_of_int %in% c("1.win5", "2.win6", "3.win9") ),
    aes(
      x=pos/1e6,
      y=-log10(rnp.clean)),
    size = 0.9,
    alpha = 0.4
    ) +
  #geom_vline(data=melt(data.frame(starts_bp, ends_bp, window_of_int = c("1.win5", "2.win6", "3.win9"))),
  #           aes(xintercept = value/1e6),
  #           color = "red",
  #           size = 0.8) +
  geom_vline(data=melt(data.frame(starts_snp, ends_snp, window_of_int = c("1.win5", "2.win6", "3.win9"))),
              aes(xintercept = value/1e6),
              color = "purple",
              size = 0.8) +
  geom_point(data = filter(outs, window_of_int %in% c("1.win5", "2.win6", "3.win9")),
             aes(x=pos/1e6,
                 y=-log10(rnp.clean),
                 fill = Annotation),
             size = 2,
             shape = 23
             ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  #geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  geom_label_repel(data = filter(outs, window_of_int %in% c("1.win5", "2.win6", "3.win9"),
                                Annotation %in% c(
                                  #"3_prime_UTR_variant",
                                  #"5_prime_UTR_variant",
                                  "missense_variant")),
                    aes(x=pos/1e6, 
                        y=-log10(rnp.clean),
                        label= Gene_Name
                        ), box.padding = 1.0, max.overlaps = Inf, size = 1.7, color = "blue") +
  theme_classic() +
  theme(legend.pos = "bottom") +
  facet_wrap(~window_of_int, scales = "free_x") -> p_val_dist

ggsave(p_val_dist, file = "p_val_dist.pdf", w= 9, h = 3)





