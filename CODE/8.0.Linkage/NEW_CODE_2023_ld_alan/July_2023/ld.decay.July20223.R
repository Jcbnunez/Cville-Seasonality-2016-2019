### LD decay

library(tidyverse)
library(magrittr)
library(reshape2)


###

ls_o = get(load("/project/berglandlab/alan/pairwise_ld_window_50000_10000.Rdata"))

# in2lt  boundaries
# 2L	2225744	13154180	2Lt
####
ls_o %>%
mutate(mid_delta = abs(mid2-mid1)) %>%
mutate(inv.stA =
case_when(
start1 < 2225744 | stop1 > 13154180 ~"noINV",
start1 >= 2225744 & stop1 <= 13154180 ~"INV"),
inv.stB =
case_when(
start2 < 2225744 | stop2 > 13154180 ~"noINV",
start2 >= 2225744 & stop2 <= 13154180 ~"INV")) %>% 
mutate(inv.ctr = paste(inv.stA, inv.stB, sep = "_")) %>%
mutate(inv.ctr = case_when(inv.ctr %in% c("INV_noINV","noINV_INV") ~ "noINV_INV",
TRUE ~ inv.ctr )  )->
ls_o.delta
ls_o.delta$inv.ctr %>% table
#######

ls_o.delta %>%
group_by(mid_delta, inv.ctr ) %>%
summarize(
mr2 = mean(meanR2, na.rm = T),
sdr2 = sd(meanR2, na.rm = T)
) %>%  
ggplot(aes(
x=(mid_delta)/1e6,
y=mr2,
ymin=mr2-sdr2,
ymax=mr2+sdr2
)) +
xlab("Mean Distance among 10kb windows") +
ylab("Mean r2 among 10kb windows") +
theme_bw() +
ylim(-0.01, 0.15) +
facet_grid(~inv.ctr) +
geom_ribbon(alpha = 0.2) +
geom_line(size = 1.2) +
geom_smooth(method = "lm", size = 1.1, color = "blue", linetype = "dashed") ->
deltar2mm
ggsave(deltar2mm, file = "deltar2mm.pdf", w = 8, h = 4)

###
ls_o.delta %>%
group_by(mid_delta, inv.ctr ) %>%
summarize(
mr2 = mean(meanR2, na.rm = T),
sdr2 = sd(meanR2, na.rm = T)
) %>%  
ggplot(aes(
x=(mid_delta)/1e6,
y=mr2,
ymin=mr2-sdr2,
ymax=mr2+sdr2,
color=inv.ctr
)) +
xlab("Mean Distance (Mb) among 10kb windows") +
ylab("Mean r2 among 10kb windows") +
theme_bw() +
ylim(-0.01, 0.14) +
geom_line(size = 1.2) +
geom_smooth(method = "lm", size = 0.9, linetype = "dashed") ->
deltar2mm_1plot
ggsave(deltar2mm_1plot, file = "deltar2mm_1plot.pdf", w = 5, h = 4)

#lms
ls_o.delta %>%
group_by(mid_delta, inv.ctr ) %>%
summarize(
mr2 = mean(meanR2, na.rm = T),
sdr2 = sd(meanR2, na.rm = T)
) %>%  
ggplot(aes(
x=(mid_delta)/1e6,
y=mr2,
ymin=mr2-sdr2,
ymax=mr2+sdr2,
color=inv.ctr
)) +
xlab("Mean Distance (Mb) among 10kb windows") +
ylab("Mean r2 among 10kb windows") +
theme_bw() +
ylim(-0.01, 0.14) +
#geom_line(size = 1.2) +
geom_smooth(method = "lm", size = 0.9, linetype = "dashed") ->
deltar2mm_1plot.lm
ggsave(deltar2mm_1plot.lm, file = "deltar2mm_1plot.lm.pdf", w = 5, h = 4)
