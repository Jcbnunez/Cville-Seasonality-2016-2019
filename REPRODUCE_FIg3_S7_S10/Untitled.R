gff.msp.300 %>%
  filter(isoform == "isoform D") %>%
  select(X4, X5) %>%
  mutate(dist= X5-X4) -> exons

exons %<>% 
  mutate(exon.of.int = case_when(X4 < 5192177 & X5 > 5192177 ~ "yea",
                                 TRUE ~ "nah"))

exons %>% data.frame()