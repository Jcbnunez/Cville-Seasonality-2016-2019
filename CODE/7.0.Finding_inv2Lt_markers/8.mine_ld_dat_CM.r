## plot LD correlation analysis
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
library(forcats)

# import master metadata
snp_guide_file <- fread("/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/Inv2L_markers_to_use_CM.txt")


#load in the ld data

ld_in_files = list()
read_in_folder <- "/project/berglandlab/Dmel_genomic_resources/Inversions/CM_2LT_markers/Rank99"

files <- system( paste("ls", read_in_folder, sep = " " ) , intern = T)


for(i in 1:length(files)){
  
  print(paste(i, "of", length(files), "or" ,i/length(files)*100, "percent", sep = " "))
  
  tmp <- fread(paste(read_in_folder, files[i], sep = "/"), sep = "\ " )
  names(tmp)[1] = c("id_tag_rank")
  
  tmp$id_tag_rank =  gsub("\\t", "", tmp$id_tag_rank)
    
  tmp %<>%
    separate(id_tag_rank, into = c("id_tag", "rank"), sep = "\\|")
  
  ld_in_files[[i]] = tmp
  
}

ld_df = do.call(rbind, ld_in_files )

#Find stragglers
#snp_guide_file_99 = snp_guide_file[which(snp_guide_file$cor_rank > 0.98),]
#data.frame(sample_ids = unique(ld_df$id_tag),
#           status= "done") %>%
#  full_join( data.frame(sample_ids = snp_guide_file_99$dm6_ID,
#                        roenumber = rownames(snp_guide_file_99) )) ->
#  sanity_check

save(ld_df, file = "ld_df.Cville.rank99.Rdata")

#summarize
ld_df %>%
  filter(R2 > 0.98)


#load in data
load("./ld_df.Cville.rank99.Rdata")

#modify rank variable
ld_df$rank = as.numeric(ld_df$rank)

## summarize
ld_df %>%
  group_by(rank) %>%
  summarise(meanR2 = mean(R2),
            sdR2 = sd(R2)
  )

in2lt_beg=2225744	
in2lt_end=13154180

ld_df %>%
  filter(R2 > 0.7) %>%
  ggplot(
    aes(
      x=as.numeric(BP_B),
      y=as.numeric(R2),
      color = as.numeric(rank)
    )
  ) + 
  geom_point() +
  geom_vline(xintercept = in2lt_beg) +
  geom_vline(xintercept = in2lt_end) +
  facet_wrap(~id_tag)->
  r2_cm

ggsave(r2_cm,
       file = "cm_r2.png",
       width = 20,
       height = 20
)

###################
################### PART 2 ---> compare and contrast

#### Finding the correlation between SNPs in CM and DGRP
load("/project/berglandlab/Dmel_genomic_resources/Inversions/DGRP_2lt_Markers/ld_df.rank90_99.Rdata")
dgrp_ld_df = ld_df
dgrp_ld_df %>%
  filter(rank > 0.98) %>%
  select(id_tag, SNP_B, R2) -> dgrp_ld_df_r
names(dgrp_ld_df_r)[3] = "R2_dgrp"


load("./ld_df.Cville.rank99.Rdata")
cm_ld_df = ld_df
cm_ld_df %>%
  select(id_tag, SNP_B, R2) -> cville_ld_df_r
names(cville_ld_df_r)[3] = "R2_cville"

#######
left_join(cville_ld_df_r, dgrp_ld_df_r) %>% 
  .[complete.cases(.),] ->
  cm_and_dgrp_r2_2lt

############3
###############################
################### PART 3 ---> the joint object

#save(cm_and_dgrp_r2_2lt,
#     file = "combined_cm_and_dgrp_r2_2lt.Rdata")
load("./combined_cm_and_dgrp_r2_2lt.Rdata")

cm_and_dgrp_r2_2lt %>%
  .[which(.$R2_cville > 0.9 & .$R2_dgrp > 0.9),] ->  
  high_ld_markers

c(high_ld_markers$id_tag, 
  high_ld_markers$SNP_B) %>%
  unique ->
  in2lt_ld_informative_markers

cm_and_dgrp_r2_2lt %>%
  filter(id_tag %in% in2lt_ld_informative_markers &
           SNP_B  %in% in2lt_ld_informative_markers) %>%
  separate(SNP_B, 
           remove = F, 
           into = c("chr", 
                    "pos", 
                    "type")) %>% 
  separate(id_tag, 
           remove = F, 
           into = c("chr_o", 
                    "pos_o", 
                    "type_o")) ->
  cm_and_dgrp_r2_2lt_ld
  
cm_and_dgrp_r2_2lt_ld %>%
  group_by(id_tag) %>%
  summarise(median = quantile(R2_cville, 0.5),
            iqr25 = quantile(R2_cville, 0.25),
            iqr75 = quantile(R2_cville, 0.75)
            ) ->
  cm_and_dgrp_r2_2lt_ld_stats

cm_and_dgrp_r2_2lt_ld_stats %>%
ggplot(aes(
  x= fct_reorder(id_tag, median),
  y= median,
  ymin = iqr25,
  ymax = iqr75
)) +
  geom_errorbar(width = 0.1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  theme(axis.text.y = element_text(size = 6)) +
  geom_point() +
  xlab("SNP ID") +
  ylab(expression(Median(r^2))) +
  coord_flip() ->
    r2_dist_per_marker

ggsave(r2_dist_per_marker,
       file = "r2_dist_per_marker.pdf",
       width = 4,
       height = 6
)

#### select the predictive snps
cm_and_dgrp_r2_2lt_ld_stats %>%
  filter(median >= 0.8) %>%
  select(id_tag) ->
  select_inv2L_markers

write.table(select_inv2L_markers, 
            file = "./in2lt_ld_informative_markers.txt", 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)  

  
## plot SNPs  
in2lt_beg=2225744	
in2lt_end=13154180

cm_and_dgrp_r2_2lt %>%
  filter(id_tag %in% select_inv2L_markers$id_tag) %>%
  filter(SNP_B %in% select_inv2L_markers$id_tag) %>%
  filter(SNP_B != id_tag) %>%
  separate(SNP_B, 
           remove = F, 
           into = c("chr", 
                    "pos", 
                    "type")) %>% 
  separate(id_tag, 
             remove = F, 
             into = c("chr_o", 
                      "pos_o", 
                      "type_o")) %>% 
 ggplot(aes(
  x= as.numeric(pos),
  y= R2_cville,
  #color = variable
)) +
  geom_point(size = 1.1, shape = 21, fill = "grey") +
  geom_vline(xintercept = in2lt_beg) +
  geom_vline(xintercept = in2lt_end) +
  geom_vline(aes(xintercept = as.numeric(pos_o)), linetype = "dashed", color = "red") +
  facet_wrap(~id_tag) +
  ylim(0,1) +
  theme_bw() ->
  inf_mark_r2_dgrp_cm

ggsave(inf_mark_r2_dgrp_cm,
       file = "inf_mark_r2_dgrp_cm.png",
       width = 12,
       height = 7
)

  




