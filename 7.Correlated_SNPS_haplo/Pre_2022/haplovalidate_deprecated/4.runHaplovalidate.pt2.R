#####
## haplovalidate
library(haploReconstruct)  
library(psych)
library(stringr)
library(data.table)
library(haplovalidate)
library(tidyverse)
library(magrittr)
library(ggrepel)

load("/project/berglandlab/summarized_dest_glm/glm.out.VA_ch_0.Rdata")

glm.out %>%
  filter(
    mod == "aveTemp",
    chr == "2L")  %>%
  select(chr, pos, b) ->
  glm
names(glm)[3] = "score"
glm$score = as.numeric(glm$score)

glm.out %>%
  filter(
    rnp.clean<0.001,
    mod == "aveTemp",
    chr == "2L")  %>%
  select(chr, pos, b) ->
  glm_out
names(glm_out)[3] = "score"
glm05$score = as.numeric(glm_out$score)
glm_out %<>%
mutate(SNP_id = paste(chr, pos))

##
load("./cands.all.rdata")
#comp_cases_id_dedupl_sort %<>%
#  .[,-which(names(comp_cases_id_dedupl_sort) == "SNP_id")]

ex_dat = comp_cases_id_dedupl_sort

ex_dat %<>% 
  mutate(SNP_id = paste(chr, pos))
  
ex_dat %<>%
  filter(SNP_id %in% glm_out$SNP_id) 

ex_dat %<>%
  .[,-which((names(ex_dat) == "SNP_id"))]

dat_filtered=initialize_SNP_time_series(chr=ex_dat$chr, 
                                        pos=ex_dat$pos,
                                        base.freq=ex_dat$basePops, 
                                        lib.freqs=ex_dat[,7:ncol(ex_dat), with=FALSE],
                                        pop.ident=rep(1:3,5), 
                                        pop.generation=c(rep(0,3),
                                                         rep(3,3),
                                                         rep(6,3),
                                                         rep(9,3),
                                                         rep(12,3)), 
                                        use.libs=c(rep(TRUE,15),rep(FALSE, 0)),
                                        minrepl = 3,
                                        minfreqchange = 0.0,
                                        max.minor.freq = 0.90,
                                        min.minor.freq = 0.05,
                                        winsize = 10000000)

dat_reconst =reconstruct_hb(dat_filtered, 
                           chrom="2L",
                           min.cl.cor = 0.70)

#save(dat_reconst,
#     file = "dat_reconst.Rdata")

#load("./dat_reconst.Rdata")

summary(dat_reconst) %>%
  .[order(-.$len.pos),] ->
  #filter(spos >= 2225744 ,
  #       epos <=  13154180) -> 
  hb_recons_table

hb_recons_table %>%
  ggplot(aes(
    ymin= id-0.2,
    ymax = id+0.2,
    xmin=spos,
    xmax=epos,
    fill = n.marker
  )) +
  geom_vline(xintercept = 2225744,
             linetype = "dashed") +
  geom_vline(xintercept = 13154180,
             linetype = "dashed") +
  geom_rect(
            color = "black") -> hapls_rec

ggsave(hapls_rec,
       file = "hapls_rec.pdf")

#cluster 3	
for(i in 1:11){
pdf(paste("win",i,"pdf", sep = "."))
plot_marker_trajectories(dat_reconst, hbr_id=i)
dev.off()
}

cluster_markers = list()
for(i in 1:11){
cluster_markers[[i]] =
  data.frame(
  cluster = i,
  chr = "2L",
  pos = markers(dat_reconst, hbr_id=i)
  )
}
cluster_markers = do.call(rbind, cluster_markers)

#### add annotations
load("../4.GLM_snps/glm_Cats.Rdata")

left_join(cluster_markers, glm_Cats[which(glm_Cats$mod == "aveTemp"),]) ->
  cluster_markers_annotated

cluster_markers_annotated %>%
  ggplot(aes(
    x=as.factor(cluster),
    y=abs(as.numeric(b)),
    color=as.factor(perm)
  )) + 
  geom_boxplot() +
  facet_wrap(.~mod) ->
  box_plot_clusters

ggsave(box_plot_clusters,
       file = "box_plot_clusters.pdf",
       width = 9,
       height = 6)

###
cluster_markers_annotated %>%
  filter(perm == 0) %>% 
  group_by(Effect, Symbol, cluster) %>%
  summarize(N = n()) %>%
  ggplot(aes(
    x=Symbol,
    fill=Effect,
    y=N
  )) +
  coord_flip() +
    geom_bar(stat = "identity") +
  facet_wrap(.~cluster, scales = "free_y")->
    bar_plot_clusters
  
  ggsave(bar_plot_clusters,
         file = "bar_plot_clusters.pdf",
         width = 12,
         height = 9)

  cluster_markers_annotated$Symbol %>% unique  %>%
    as.data.frame()
  
  
##3 bring in LD data
##
load("/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/ld_tempOut_plus_inv.Cville.rank99.Rdata")

ld_df %>%
  filter(CHR_A == "2L") %>% 
  filter(BP_A %in% cluster_markers_annotated$pos) %>% 
  mutate(delta_bp = abs(BP_A-BP_B)) %>%
  filter(delta_bp > 0) %>% 
  ggplot(
    aes(
      x=delta_bp,
      y=R2
    )
  ) + geom_smooth() ->
  all_SNPS_2l_r2

ggsave(all_SNPS_2l_r2,
       file = "all_SNPS_2l_r2.png",
       width = 6,
       height = 6)


ld_df %>%
  filter(CHR_A == "2L") %>%
  filter(BP_A %in% cluster_markers_annotated$pos) %>% 
  filter(BP_B %in% cluster_markers_annotated$pos) %>% 
  mutate(pos = BP_A) %>% 
  left_join(cluster_markers_annotated[which(cluster_markers_annotated$perm == 0),] ) %>% 
  ggplot(aes(
    x=(abs(BP_A-BP_B)),
    y=R2,
    color = as.factor(cluster)
  )) + 
  scale_x_continuous(trans='log10') +
    geom_point(size =1) ->
    points_R2_clusters
  
  ggsave(points_R2_clusters,
         file = "points_R2_clusters.pdf",
         width = 6,
         height = 6)
  
### R> 0.75
  cluster_markers_annotated[which(cluster_markers_annotated$perm == 0),] -> annot_filt

  annot_filt %>%
    select(pos, Symbol) -> base_annot
  
  A_annot = base_annot
  names(A_annot)[1:2] = c("BP_A", "A_annot")
  B_annot = base_annot
  names(B_annot)[1:2] = c("BP_B", "B_annot")
  
  ld_df %>%
    filter(CHR_A == "2L") %>%
    filter(BP_A %in% cluster_markers_annotated$pos) %>% 
    filter(BP_B %in% cluster_markers_annotated$pos) %>% 
    mutate(pos = BP_A) %>% 
    left_join( annot_filt) -> markers_of_interest
 
  markers_of_interest   %>% 
    filter(R2 >= 0.70,
           cluster == 1) %>%
    left_join(A_annot) %>% 
    left_join(B_annot) %>%
    mutate(delta_bp = abs(BP_A-BP_B)) %>%
    filter(delta_bp > 0) %>%
    mutate(ld_pair = paste(A_annot, B_annot, sep = ":")) %>%
    group_by(ld_pair) %>%
    slice_head() %>%
    ggplot(aes(
      x=delta_bp,
      y=R2,
      color = as.factor(cluster),
      label = ld_pair
    )) + 
    scale_x_continuous(trans='log10') +
    geom_point(size =3) +
    geom_text_repel() ->
    points_R2_cluster1
  
  ggsave(points_R2_cluster1,
         file = "points_R2_cluster1.pdf",
         width = 12,
         height = 6)
  
  
  
  markers_of_interest   %>% 
    filter(Symbol == "poe") %>%
    left_join(A_annot) %>% 
    left_join(B_annot) %>%
    mutate(delta_bp = abs(BP_A-BP_B)) %>%
    filter(delta_bp > 0) %>%
    mutate(ld_pair = paste(A_annot, B_annot, sep = ":")) %>%
    group_by(ld_pair) %>%
    slice_head() %>%
    ggplot(aes(
      x=delta_bp,
      y=R2,
      color = as.factor(cluster),
      label = ld_pair
    )) + 
    scale_x_continuous(trans='log10') +
    geom_point(size =3) +
    geom_text_repel() ->
    points_R2_poe
  
  ggsave(points_R2_poe,
         file = "points_R2_poe.pdf",
         width = 12,
         height = 6)
  
  