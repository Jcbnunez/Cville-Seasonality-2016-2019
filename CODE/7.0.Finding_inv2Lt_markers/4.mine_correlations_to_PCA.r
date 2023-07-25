## find the SNPs woth perfect correlation to inversion status
## 

#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)


load("./PCA_obj_2LT_dimdesc.Rdata")


PCA_obj_2LT_dimdesc$Dim.1$quanti %>%
  .[complete.cases(.),] %>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, into = c("chr","pos","Type"), remove = F) %>% 
  mutate(cor_rank = round( abs(correlation), ),
         p.value.bonfe = p.adjust(p.value, method = "fdr")) ->
  correlation_info

correlation_info %>%
  group_by(cor_rank) %>%
  summarise(N = n())

in2lt_beg=2225744	
in2lt_end=13154180

correlation_info %>%
  ggplot(
    aes(
      x=as.numeric(pos),
      y=-log10(p.value.bonfe),
      color=cor_rank
    )
  ) + 
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  xlab("2L Position (bp)") +
  ylab( expression(-log[10](italic(P)[Bonferroni])))  +
  geom_vline(xintercept = in2lt_beg) +
  geom_vline(xintercept = in2lt_end) ->
  p_val_dist

ggsave(p_val_dist,
       file = "p_val_dist.png", 
       width = 6,
       height  = 3
       )

#Generate final file of hits

correlation_info %>%
  .[which(.$p.value.bonfe < 0.01 &
          abs(.$correlation) >= 0.5),] %>%
  group_by(cor_rank) %>%
  summarize(N= n() ) %>%
  ggplot(
    aes(
      x=cor_rank,
      y=N
    )
  ) +
  geom_point() +
  ylab("Number of SNPs") +
  xlab(expression(absolute(rho))) ->
  cor_v_N

ggsave(cor_v_N,
       file = "cor_v_N.pdf", 
       width = 4,
       height  = 4
)

###Print data

correlation_info %>%
filter(Type == "SNP") %>%
  .[which(.$p.value.bonfe < 0.01 &
            abs(.$correlation) >= 0.5),] ->
  inv2Lt_markers

write.table(inv2Lt_markers, 
            file = "./inv2L_informative_markers_Dm3.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = TRUE)

