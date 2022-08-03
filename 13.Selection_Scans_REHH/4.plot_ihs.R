### libraries
rm(list = ls())

library(data.table)
library(SeqArray)
library(rehh)
library(tidyverse)
library(foreach)
library(car)
library(DescTools)
library(doMC)
registerDoMC(2) ## using an i-job ##ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

final.windows.pos = 
  data.frame(win.name = c("win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6" ),
             mid = c(3.0, 4.67, 5.15, 6.2, 6.8 , 9.6)
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


setwd("/scratch/yey2sn/Overwintering_ms/13.Selection_REHH")

load("CM.IHS.HH.INV.out.Rdata")
wgscan.ihs$ihs %>% 
  as.data.frame() ->
  ihs.dat

setDT(ihs.dat)
###############
### windows ###
###############
# generate a master index for window analysis
### define windows
win.bp <- 1e5
step.bp <- 5e4

setkey(ihs.dat, "CHR")

## prepare windows
wins <- foreach(chr.i=c("2L"
                        #,"2R", "3L", "3R"
                        ),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- ihs.dat %>%
                    filter(CHR == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$POSITION), to=max(tmp$POSITION)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$POSITION), to=max(tmp$POSITION)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]

dim(wins)
#####
#####

ihs.win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- ihs.dat %>%
    filter(CHR ==  wins$chr[win.i]) %>%
    filter(POSITION > wins$start[win.i] & POSITION <  wins$end[win.i])
  
  data.frame(
    CHR =  wins$chr[win.i],
    MIDPOS = c( wins$start[win.i] + wins$end[win.i])/2,
    MeanIHS.abs = mean(abs(win.tmp$IHS), na.rm = T),
    MeanIHS.raw = mean((win.tmp$IHS), na.rm = T),
    MeanP.val = mean((win.tmp$LOGPVALUE), na.rm = T)
  )
  
}

ihs.win.out = win.out

mean.abs.ihs = mean(abs(ihs.dat$IHS), na.rm = T )

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = 6671911, col = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 6202053, col = "blue", linetype = "dashed") +
  geom_vline(xintercept = 2225744) +
  geom_vline(xintercept = 13154180) +
 geom_rect(data = final.windows.pos,
           aes(xmin=start, xmax = end,
               ymin = -0.6, ymax = 0.6), 
           alpha = 0.7, fill = "gold") +
  geom_line(data = ihs.win.out,
            aes(
              x=MIDPOS,
              y=MeanIHS.raw
            ))  +
  ylab("Mean iHS (per 0.1 Mb window)") +
  ggtitle("Inverted Karyotypes") +
  xlab("Genomic Position (Mb)") ->
  ihs_line.inv

ggsave(ihs_line.inv, file = "ihs_line.inv.pdf", w = 9, h = 3)
ggsave(ihs_line.inv, file = "ihs_line.inv.png", w = 9, h = 3)


#####

load("rehh.INV.object.CM.Rdata")

hh@positions[ which(names(hh@positions) %in%  c("2L_6671911_SNP", "2L_6202053_SNP") )]

furcation <- calc_furcation(hh,
                            mrk = "2L_6202053_SNP")

pdf("reeh.pdf")
plot(furcation)
dev.off()


