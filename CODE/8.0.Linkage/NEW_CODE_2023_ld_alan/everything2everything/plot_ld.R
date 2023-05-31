## Beging by pasting file

#.  system("cp /project/berglandlab/alan/pairwise_ld_window_50000_10000.Rdata ./")

### libraries
  library(ggplot2)
  library(data.table)
  library(viridis)
  library(tidyverse)
  library(magrittr)
  library(gmodels)
  
### load data
  load("./pairwise_ld_window_50000_10000.Rdata")

#### Quantify data -->
  final.windows.pos1 = 
    data.frame(win.name = 
    c("bkpt", "win", "win", "win", "win", "win", "win", "bkpt" ),
               mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
               chr = "2L") %>%
    mutate(start1 = (mid-0.2)*1e6,
           stop1  = (mid+0.2)*1e6)
  setDT(final.windows.pos1)
  ##
  final.windows.pos2 = 
    data.frame(win.name = 
                 c("bkpt", "win", "win", "win", "win", "win", "win", "bkpt" ),
               mid = c(2.2, 3.0, 4.67, 5.12, 6.2, 6.8 , 9.6, 13.1),
               chr = "2L") %>%
    mutate(start2 = (mid-0.2)*1e6,
           stop2 = (mid+0.2)*1e6)
  setDT(final.windows.pos2)

####
####
setkey(final.windows.pos1, start1, stop1)
foverlaps(o, final.windows.pos1) %>% .$win.name -> win1.name

setkey(final.windows.pos2, start2, stop2)
foverlaps(o, final.windows.pos2) %>% .$win.name -> win2.name

o %>%
  mutate(win1name = win1.name) %>%
  mutate(win2name = win2.name) %>%
  mutate(joint_nam = paste(pmin(win1name,win2name), 
                           pmax(win1name,win2name), 
                           sep = ";"))->
  o.mapped

o.mapped %<>%
  mutate(outside = case_when(mid1 > 13154180 &  mid2 > 13154180 ~ "Out"))

o.mapped %>%
  group_by(joint_nam) %>%
  summarize(meanR2 = mean(meanR2, na.rm = T),
            sdR2 = sd(meanR2, na.rm = T),
            ) %>% as.data.frame()
o.mapped %>%
  group_by(outside) %>%
  summarize(meanR2 = mean(meanR2, na.rm = T),
            sdR2 = sd(meanR2, na.rm = T),
  ) %>% as.data.frame()
###  
### pad empty spaces
  grid <- data.table(expand.grid(1:max(c(o$win1, o$win2)), 1:max(c(o$win1, o$win2))))
  setnames(grid, names(grid), c("win1", "win2"))

  grid <- grid[win1<win2]
  setkey(o, win1, win2)
  setkey(grid, win1, win2)

  o2 <- merge(o[poolOnly==T][rnp.thr==1], grid, all=T)
  o2[is.na(meanR2), meanR2:=-.01]
  
  table(o$win1>o$win2)
  table(o2$win1>o2$win2)

### plot
  p1 <- ggplot(data=o2[abs(win1-win2)>10][meanR2>0], 
               aes(x=mid1/1e6, y=mid2/1e6, fill=meanR2)) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = 20, linetype = "solid") +
    geom_hline(yintercept = 2225744/1e6, linetype = "solid") +
    geom_hline(yintercept = 13154180/1e6, linetype = "solid") +
    geom_vline(xintercept = 4.7, linetype = "solid") +
    geom_vline(xintercept = 5.2, linetype = "solid") +
    geom_vline(xintercept = 6.1, linetype = "solid") +
    geom_vline(xintercept = 6.8, linetype = "solid") +
    geom_vline(xintercept = 9.6, linetype = "solid") +
    #geom_raster(data=o2[win1!=win2][meanR2<0], aes(x=mid1, y=mid2), fill="grey39", alpha=.85) +
                geom_raster() + 
    scale_fill_stepsn(#colours=c("steelblue","springgreen3","orange2","firebrick"),
                      colours=rev(rainbow(5)),
                      breaks = seq(from = 0, to = 0.10, by = 0.025),
                      limits=c(0,0.10)) +
    xlim(0,20) +
    ylim(0,20) +
    #scale_fill_gradient2(midpoint = 0.0005, 
    #                     low = "steelblue", high = "firebrick") +
        theme_classic() 

  ggsave(p1, file="./ld_mat_r.pdf", h=12, w=12)


#ggplot(data=o[win1!=win2], aes(x=win1, y=win2, fill=meanR2)) + geom_tile() + scale_fill_viridis(option="F") + #facet_grid(poolOnly~rnp.thr)
#