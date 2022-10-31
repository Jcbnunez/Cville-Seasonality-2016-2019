### Code to process ld data
### 

library(tidyverse)
library(data.table)
library(magrittr)
library(forcats)

###
###
system(intern = T, "ls ld_plink | grep -v 'log' | grep -v 'nosex'") -> ld_files
ld_list = list()

for(i in 1:length(ld_files)){

print(i)
  
tmp <- fread(paste("./ld_plink/", ld_files[i], sep = ""))

mins <- apply(tmp[,c("BP_A","BP_B")], 1, min)
maxs <- apply(tmp[,c("BP_A","BP_B")], 1, max)

tmp %<>%
  mutate(pair_id = paste(mins,maxs, sep = "_" ) )

ld_list[[i]] = tmp
}

ld_df = do.call(rbind, ld_list)


##########


system(intern = T, "ls corrOut") -> cor_files
cor_list = list()

for(i in 1:length(cor_files)){
  
  print(i)
  
load(paste("./corrOut/", cor_files[i], sep = ""))
tmp2 <- out_df

mins <- apply(tmp2[,c("pos1","pos2")], 1, min)
maxs <- apply(tmp2[,c("pos1","pos2")], 1, max)

tmp2 %<>%
  mutate(pair_id = paste(mins,maxs, sep = "_" ) )

cor_list[[i]] = tmp2
}

cor_df = do.call(rbind, cor_list)

###
###
joint_df = left_join(ld_df, cor_df, by = "pair_id" )

joint_df %>% 
  .[complete.cases(.),] ->
  joint_df_cc

save(joint_df_cc, file = "joint_df_cc.Rdata")
######
######

load("./joint_df_cc.Rdata")
joint_df_cc %>% head

selected_snps <-
unique(unique(joint_df_cc$SNP_A),
unique(joint_df_cc$SNP_B))

selected_snps = gsub("_SNP", "", selected_snps)

### load annotations
load("/project/berglandlab/DEST_Charlottesville_TYS/Annotation_Dmel6_SNAPE.Rdata")
annotation_dmel_sp %>% head

annotation_dmel_sp %>%
  filter(SNP_id %in% selected_snps) %>%
  select(chr, pos, Symbol, FBgn) %>%
  mutate(SNP_A = paste(chr,pos, "SNP", sep = "_")) %>% 
  select(Symbol, FBgn, SNP_A)->
  annot_A
names(c) = c("Symbol_A", "FBgn_A", "SNP_A" )

annotation_dmel_sp %>%
  filter(SNP_id %in% selected_snps) %>%
  select(chr, pos, Symbol, FBgn) %>%
  mutate(SNP_B = paste(chr,pos, "SNP", sep = "_")) %>% 
  select(Symbol, FBgn, SNP_B)->
  annot_B
names(annot_B) = c("Symbol_B", "FBgn_B", "SNP_B" )

####
left_join(joint_df_cc, annot_A) %>%
  left_join(annot_B) ->
  joint_df_cc_annot

####
joint_df_cc_annot$pair_id %>% table %>% table

joint_df_cc_annot %>% 
  group_by(pair_id) %>%
  slice_head() ->
  joint_df_cc_annot_dedup

joint_df_cc_annot_dedup$pair_id %>% table %>% table

#####
save(joint_df_cc_annot_dedup, file = "joint_df_cc_annot_dedup.Rdata")
######
######
load("./joint_df_cc_annot_dedup.Rdata")

joint_df_cc_annot_dedup %>%
  as.data.frame %>% 
  head

joint_df_cc_annot_dedup$pos1 = as.numeric(joint_df_cc_annot_dedup$pos1)
joint_df_cc_annot_dedup$pos2 = as.numeric(joint_df_cc_annot_dedup$pos2)

apply(joint_df_cc_annot_dedup[,c("pos1", "pos2")], 1, min ) -> mins
apply(joint_df_cc_annot_dedup[,c("pos1", "pos2")], 1, max ) -> maxs

joint_df_cc_annot_dedup %<>%
  as.data.frame() %>%
  mutate(Pos_min = as.numeric(mins),
         Pos_max = as.numeric(maxs)) 
 
#### Plots

joint_df_cc_annot_dedup %>% 
  ggplot(aes(R2)) +
  geom_histogram() ->
  R2_hist

joint_df_cc_annot_dedup %>% 
  ggplot(aes(correlation)) +
  geom_histogram() ->
  correlation_hist

ggsave(correlation_hist,
       file = "correlation_hist.pdf",
       width = 5,
       height = 4)


## correlation plot
joint_df_cc_annot_dedup %>%
  ggplot(aes(
    y=R2,
    x=abs(correlation)
  )) +
  geom_point() +
  geom_smooth(method = "lm")->
  cor_v_r2

ggsave(cor_v_r2,
       file = "cor_v_r2.png",
       width = 5,
       height = 4)


cor.test(joint_df_cc_annot_dedup$R2, abs(joint_df_cc_annot_dedup$correlation))


### threshold plot

joint_df_cc_annot_dedup %>% 
  filter(abs(correlation) > 0.6,
         R2 > 0.6,
         Pos_min < 1e7,
         Pos_max < 1e7) %>% 
  ggplot(aes(
    x=(Pos_min/1e6),
    y=(Pos_max/1e6),
    color = R2,
    size = R2
  )) +
  geom_abline(intercept = 0, slope = 1, size = 0.7, linetype = "dashed") +
  geom_jitter(shape = 20, 
             alpha = 0.7) +
  theme_classic() +
  scale_size(range = c(0.05, 5)) +
  ggtitle("R2 for SNP pairs with PoolSeq", subtitle = "R2 > 0.6 and |Cor| > 0.6") +
  xlab("SNP A (Mb)") +
  ylab("SNP B (Mb)") +
  scale_color_gradient(low = "blue", high = "red", na.value = NA) ->
  R2_point_plot

ggsave(R2_point_plot,
       file = "R20.5_point_plot.pdf",
       width = 5,
       height = 4)

#####
joint_df_cc_annot_dedup %>% 
  filter(abs(correlation) > 0.6,
         R2 > 0.6,
         bp_dist > 200000,
         bp_dist < 1000000) -> long_distance_markers

c(long_distance_markers$Symbol_A, long_distance_markers$Symbol_B) %>% 
  table %>%
  as.data.frame() %>%
  filter(Freq > 5) %>%
  ggplot(aes(y=fct_reorder(`.`, Freq),
         x=Freq)) +
    theme(text = element_text(size=10)) +
    geom_bar(stat = "identity") ->
    genes_involved
  
  ggsave(genes_involved,
         file = "genes_involved.pdf",
         width = 5,
         height = 4)
  
#### intertaction table
  long_distance_markers %>%
    mutate(gene_combo = paste(Symbol_A, Symbol_B, sep = "|")) %>%
    .$gene_combo %>%
    table %>% 
    as.data.frame() -> 
    table_of_interactions

  table_of_interactions[grep("-", table_of_interactions$`.`, invert = T),]  ->
    table_of_interactions
  
  table_of_interactions[grep("smal", table_of_interactions$`.`),]
  table_of_interactions[grep("stai", table_of_interactions$`.`),]
  table_of_interactions[grep("Pez", table_of_interactions$`.`),]
  
##individual interactions
  joint_df_cc_annot_dedup %>%
    filter(Symbol_A == "smal" | Symbol_B == "smal" ) %>% 
    filter(abs(correlation) > 0.6,
           R2 > 0.6) %>%
    select(Pos_min, R2, correlation) %>%
    melt(id = "Pos_min") %>%
    mutate(focal = "smal") -> smal_df
    
  joint_df_cc_annot_dedup %>%
    filter(Symbol_A == "stai" | Symbol_B == "stai" ) %>% 
    filter(abs(correlation) > 0.6,
           R2 > 0.6) %>%
    select(Pos_min, R2, correlation) %>%
    melt(id = "Pos_min") %>%
    mutate(focal = "stai") -> stai_df

### 
  rbind(smal_df, 
  stai_df) %>%
    ggplot(aes(x= Pos_min,
               y= value,
               color = focal)) +
    geom_vline(xintercept = 6179133) +
    geom_point() + 
    facet_wrap(~variable, ncol = 1)->
    genes_Cases
    
  ggsave(genes_Cases,
         file = "genes_Cases.pdf",
         width = 5,
         height = 4)
  
    