library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)


dat_in = load_in_df


### plot p
glm.out_p1$SNP_id

begin=2225744	
end=13154180

bk_output %>%
  select(!t.test.mu) %>%
  .[complete.cases(.),] %>% 
  filter(test_type == "inv>inv" |
           test_type %in% c("inv>tmp") #& affected_snp %in% glm.out_p1$SNP_id
  ) %>%
  separate(affected_snp,
           remove = F,
           into = c("chr","pos","feature"), sep = "_" ) ->
  dat_3

dat_3 %>%
  ggplot(
    aes(
      x=as.numeric(pos),
      y=log2(prediction_robustness/0.05),
      #shape = test_type,
      color = Real_R
      #color = focal_snp
    )) + 
  geom_vline(xintercept = begin) +
  geom_vline(xintercept = end) +
  geom_jitter(size = 1.5) +
  ylab(expression(log[2](P[ind]/alpha))) +
  scale_shape_manual(values = 15:16) +
  scale_color_gradientn(name = expression(R[obs]^2), colours = terrain.colors(10)) +
  facet_grid(~test_type) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") ->
  test3

ggsave(test3,
       file = "test3.png",
       width = 8,
       height = 2)

dat_3 %>%
  ggplot(
    aes(
      x= Real_R,
      y= log2(prediction_robustness/0.05),
      color = Real_R
    )) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab(expression(log[2](P[ind]/alpha))) +
  xlab(expression(R[real]^2)) +
  scale_color_gradientn(name = expression(R[obs]^2), colours = terrain.colors(10)) ->
  test3_plot_points

ggsave(test3_plot_points,
       file = "test3_plot_points.png",
       width = 5,
       height = 4)

####
#####3
##bk_output %>%
##  select(!t.test.mu) %>%
##  .[complete.cases(.),] %>% 
##  filter(test_type %in% c("tmp>tmp")) %>%
##  separate(affected_snp,
##           remove = F,
##           into = c("chr","pos","feature"), sep = "_" ) ->
##  dat_4
##
##### ComplexHeatmapComplexHeatmap ComplexHeatmapComplexHeatmapComplexHeatmapComplexHeatmap
##library(ComplexHeatmap)
##pdf("test4.pdf")
##dat_4 %>%
##  mutate(test_stat = prediction_robustness/0.05) %>%
##  dcast(affected_snp~focal_snp, value.var = "test_stat") %>% 
##  .[,-1] %>% 
##  as.matrix() %>% 
##  densityHeatmap() 
##dev.off()
##
##
##
##
##dat_4 %>%
##  ggplot(
##    aes(
##      x=as.numeric(pos),
##      y=log2(prediction_robustness/0.05),
##      color = test_type
##      #color = focal_snp
##    )) + 
##  geom_smooth() +
##  geom_point() +
##  ylab(expression(log[2](P[ind]/alpha))) +
##  geom_vline(xintercept = begin) +
##  geom_vline(xintercept = end) +
##  geom_hline(yintercept = 0, linetype = "dashed", color = "red") ->
##  test4
##
##ggsave(test4,
##       file = "test4.png",
##       width = 6,
##       height = 2)
##