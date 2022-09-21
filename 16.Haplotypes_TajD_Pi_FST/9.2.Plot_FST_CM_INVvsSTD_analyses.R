### FST analysis CM
### 
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)

gff.msp.300 = vroom("gff.msp.300.dat.txt", col_names = F, delim = "\t")
for(i in 1:dim(gff.msp.300)[1] ){
  tmp <- gff.msp.300[i,]
  
  tmp.str = str_split(tmp, pattern = "\\|")
  
  isoform = gsub('.+(isoform .).+','\\1',tmp.str[[9]][3])
  
  gff.msp.300$isoform[i] = isoform
  
}
gff.msp.300$isoform %>% unique()

gff.msp.300 %<>% 
  mutate(rank = case_when(isoform == "isoform B" ~ -0.1,
                          isoform == "isoform D" ~ -0.2,
                          isoform == "isoform E" ~ -0.3,
                          isoform == "isoform F" ~ -0.4,
                          isoform == "isoform G" ~ -0.5,
                          isoform == "isoform H" ~ -0.6,
                          isoform == "isoform I" ~ -0.7,
                          isoform == "isoform J" ~ -0.8,
                          isoform == "isoform K" ~ -0.9,
                          isoform == "isoform L" ~ -1.0,
                          isoform == "isoform M" ~ -1.1
  ))

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6,
         CHROM = "2L")

### All sites
fst.cm <- vroom("inv.vs.std.cm.fst.windowed.weir.fst", na = c("", "NA", "-nan"))

fst.cm %>%
  filter(CHROM == "2L") ->
  fst.all.2L

ggplot() +
  geom_rect(data = final.windows.pos,
            aes(xmin=start, xmax = end,
                ymin = -0, ymax = 0.65), 
            alpha = 0.7, fill = "gold") +
  
geom_line(data =fst.all.2L,  aes(
    x=c(BIN_START+BIN_END)/2,
    y=MEAN_FST
  )) +
  facet_wrap(CHROM~., scales = "free_x", ncol = 1) +
  theme_bw() +
  geom_vline(xintercept = 5170001, color = "red") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) ->
  fst.inv.std
ggsave(fst.inv.std, file ="fst.inv.std.pdf", w = 7, h = 3)

#### zoom to MSP 300
#### 

ggplot() +
  geom_rect(data = gff.msp.300,
            aes(xmin=X4, xmax = X5,
                ymin = rank/2-0.01, ymax = rank/2+0.01), 
            alpha = 0.7, fill = "red") +
  geom_text(
    data = slice_head(group_by(gff.msp.300, isoform)),
    aes(x=5.07e6, y = rank/2,
        label = isoform
    ), size = 2.5 ) +
  geom_line(data=filter(fst.all.2L, BIN_START > 5.05e6-1e5, BIN_END < 5.25e6+1e5),
            aes(
              x=(BIN_START+BIN_END)/2,
              y=MEAN_FST,
            )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  ggtitle("MSP300") +
  ylab(expression(F[ST]))  ->
  fst.inv.std.w5.1

ggsave(fst.inv.std.w5.1, file ="fst.inv.std.w5.1.pdf", w = 9, h = 6)


  