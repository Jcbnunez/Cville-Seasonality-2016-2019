library(tidyverse)
library(magrittr)
library(ggbeeswarm)
library(reshape2)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions

#######

files_to_read <- system("ls ./p_val_extractor_VA_only_granular | grep '.Rdata' ", intern = T)


collect_counts = list()
collect_enrrich = list()
for(i in 1:length(files_to_read)){
  print(files_to_read[i])
  
  load( paste("/scratch/yey2sn/Overwintering_ms/4.GML_plots/p_val_extractor_VA_only_granular/",
              files_to_read[i],
              sep = ""))
  
  collect_counts[[i]] = count_df
  collect_enrrich[[i]] = enrrich_df
  
}

counts_df = do.call(rbind, collect_counts)
counts_df %<>%
  mutate(proportion = Tresh/N)

counts_df$analysis_type %>% table

#### Now extract only the counts
counts_df %>%
  filter(analysis_type == "Annotation_enrichment") ->
  enrrich_count

dat_in=enrrich_count %>% filter(p_tresh < 0.1)
ggplot() +
  geom_violin(data=dat_in[which(dat_in$category == "Per"),],
              aes(x=p_tresh,
                  y=proportion/as.numeric(p_tresh),
                  fill = inversion_pos),
              alpha = 0.7
  ) +
  geom_point(data=dat_in[which(dat_in$category == "Obs"),],
             aes(x=p_tresh,
                 y=proportion/as.numeric(p_tresh),
                 fill = inversion_pos),
             size = 1.3,
             shape = 23,
             color = "black",
             position = position_dodge2(w = 0.95)) +
  facet_grid(Tidy_Consequence~chr, scales = "free_x", space = "free") +
  ylab(expression( paste("(%", italic(P[LRT]<alpha),")","/", alpha ) ) ) +
  xlab(expression( alpha ) )  +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "top")-> 
  p_tresh_plot

ggsave(p_tresh_plot, 
       file = "p_tresh_plot.VA.png",
       width = 12,
       height = 6)

#### ECDF

chrs = unique(enrrich_count$chr)
p_treshs = unique(enrrich_count$p_tresh)
cats = unique(enrrich_count$Tidy_Consequence)
invs = unique(enrrich_count$inversion_pos)

inner_t=list()
inner_k=list()
inner_i=list()
inner_j=list()

for(t in 1:length(invs) ) {
for(k in 1:length(cats) ) {
for(i in 1:length(p_treshs) ) {
for(j in 1:length(chrs) ) {
    
    print(paste( round(i/length(p_treshs)*100)
                ,round(j/length(chrs)*100)
                ,round(k/length(cats)*100)
                ,round(t/length(invs)*100)
                , sep = " | "))
    
  Ho <- enrrich_count$proportion[which(  enrrich_count$p_tresh == p_treshs[i] & 
                                           enrrich_count$chr == chrs[j] & 
                                           enrrich_count$Tidy_Consequence == cats[k] & 
                                           enrrich_count$inversion_pos == invs[t] & 
                                           enrrich_count$category == "Obs")] 
    
  DatProp <- enrrich_count$proportion[which(  enrrich_count$p_tresh == p_treshs[i] & 
                                              enrrich_count$chr == chrs[j] & 
                                              enrrich_count$Tidy_Consequence == cats[k] &
                                              enrrich_count$inversion_pos == invs[t] & 
                                              enrrich_count$category == "Per")] 
    
    
    #if(length(DatProp)>0){ ### protection against null values
      ecdf(DatProp) -> tmp_ecdf
      
      out_emp_p = data.frame(
        p.tresh = p_treshs[i],
        chr = chrs[j],
        Cat = cats[k],
        Inv = invs[t],
        ECDF.q = tmp_ecdf(Ho),
        t.stat = t.test(DatProp, mu = Ho)$statistic,
        t.p = t.test(DatProp, mu = Ho)$p.value
      )
      
      inner_j[[j]] = out_emp_p
    #} ### protection against null values


} #close j
    inner_i[[i]] = do.call(rbind, inner_j)
} #close i
    inner_k[[k]] = do.call(rbind, inner_i)
} #close k
    inner_t[[t]] = do.call(rbind, inner_k)
} #close t

final_ecdf_df = do.call(rbind, inner_t)
save(final_ecdf_df, file = "final_ecdf_df.Virginia.Rdata")

final_ecdf_df %>%
  ggplot(aes(
         x=(as.numeric(p.tresh)),
         y=-t.stat ,
         color=Inv)) +
  geom_hline(yintercept = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),) +
  geom_smooth(size = 0.8) +
  geom_point(size = 1.1, alpha = 0.3) +
  theme(axis.text = element_text(size = 7, angle=45)) +
  ylab(expression(paste("T test stat" )) ) +
  xlab(expression(log[10](alpha))) +
  facet_grid(chr~Cat) ->
  t.test_plot

ggsave(t.test_plot, 
       file = "t.test_plot.pdf",
       width = 12,
       height = 6)

final_ecdf_df %>%
  ggplot(aes(
    x=(as.numeric(p.tresh)),
    y=ECDF.q ,
    color=Inv)) +
  geom_hline(yintercept = 0) +
  #cale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)),) +
  geom_smooth(size = 0.8) +
  geom_point(size = 1.1, alpha = 0.3) +
  facet_grid(chr~Cat) ->
  ECDF.q.plot

ggsave(ECDF.q.plot, 
       file = "ECDF.q.plot.pdf",
       width = 12,
       height = 6)

##### Figure for main paper
##### 

final_ecdf_df %>%
  filter(p.tresh == 0.05) %>%
  mutate(Empical.p = 1-ECDF.q) %>%
  filter(Empical.p <= 0.05) %>%
  dplyr::select(chr,Cat,Inv, Empical.p ) %>%
  mutate(round_p = round(Empical.p, 3))


dat_in=enrrich_count %>% filter(p_tresh == 0.05)
ggplot() +
  geom_violin(data=dat_in[which(dat_in$category == "Per"),],
              aes(x=chr,
                  y=proportion,
                  fill = inversion_pos),
              alpha = 0.7
  ) +
  geom_point(data=dat_in[which(dat_in$category == "Obs"),],
             aes(x=chr,
                 y=proportion,
                 fill = inversion_pos),
             size = 2,
             shape = 23,
             color = "black",
             position = position_dodge2(w = 0.95)) +
  facet_grid(~Tidy_Consequence, scales = "free_x", space = "free") +
  ylab(expression( paste("%", italic(P[LRT]<alpha) ) ) ) +
  xlab(expression( alpha ) )  +
  theme_bw() +
  theme(legend.position = "top")-> 
  p_tresh_plot_p5

ggsave(p_tresh_plot_p5, 
       file = "p_tresh_plot_p5.pdf",
       width = 9,
       height = 2.8)
