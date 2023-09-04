#### plot the fst means

library(tidyverse)
library(data.table)
library(vcfR)
library(foreach)
library(magrittr)
library(SeqArray)

####

invr <- get(load("r2.invsnp.CM.Rdata"))
setDT(invr)

invr %<>%
separate(SNPid, remove = F, into = c("chr","pos"), sep = "_")
### summarize

 # generate a master index for window analysis
  ### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  
  setkey(invr, "chr")
  invr$pos = as.numeric(invr$pos)
  obj = invr ### < ---- object comes in here
  setDT(obj)
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L"),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- obj %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  
  wins[,i:=1:dim(wins)[1]]
  
  dim(wins)
  ##### 
  ##### 
  ##### 
  ##### 
  
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind"
  )%do%{
    
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.tmp <- obj %>%
    filter(chr == wins[win.i]$chr ) %>%
    filter( pos >= wins[win.i]$start & 
    		pos <= wins[win.i]$end)
    
    
    	
    win.tmp %>%
    summarize(
    chr = wins[win.i]$chr,
    pos = median(pos),
    r2.m = mean(r2, na.rm = T)) %>%
    mutate(r2.99n = sum(win.tmp$r2 > 0.99),
   		   r2.70n = sum(win.tmp$r2 > 0.70),
    	   nloc = length(win.tmp$r2)
    ) ->
    o
    
    return(o)
    
    }

#####
final.windows.pos = 
  data.frame(win.name = c(#"left", 
                          "win_3.1", "win_4.7", "win_5.1", "win_6.1", "win_6.8", "win_9.6"
                          #, "right" 
                          ),
             mid = c(#2.2, 
                     3.1, 4.7, 5.2, 6.1, 6.8 , 9.6
                     #, 13.1
                     ),
             chr = "2L"
  ) %>%
  mutate(start = (mid-0.2)*1e6 ,
         end  = (mid+0.2)*1e6  )


#####
invr %>%
mutate(Win = case_when(
    pos > 2000000 & pos < 2400000 ~ "left",
    pos > 12900000 & pos < 13300000 ~ "right",
    pos > 2800000 & pos < 3200000 ~ "w3.1",
    pos > 4470000 & pos < 4870000 ~ "w4.7",
    pos > 4920000 & pos < 5320000 ~ "w5.1",
    pos > 6000000 & pos < 6400000 ~ "w6.1",
    pos > 6600000 & pos < 7000000 ~ "w6.8",
    pos > 9400000 & pos < 9800000 ~ "w9.6")) %>%
    group_by(Win) %>%
    summarize(
    mean.r = mean(r2),
    r99c = sum(r2 > 0.99),
    r70c = sum(r2 > 0.70),
    )

#> win.out$r2.70n %>% mean
#[1] 12.06397

#####
ggplot() +
  geom_rect(data = final.windows.pos, 
  aes(xmin=start/1e6 , xmax =end/1e6, 
  ymin = 0, ymax = 0.20),
              fill = "gold", alpha = 0.6) +
geom_vline(xintercept = 2225744/1e6, color = "red") +
geom_vline(xintercept = 13154180/1e6, color = "red") +
geom_smooth(
data = win.out,
aes(
x=pos/1e6,
y=r2.m)) +
geom_line(
data = win.out,
aes(
x=pos/1e6,
y=r2.m
)) + 
ylab("mean r2 to Inv") +
xlab("Chr 2L Mbp") +
ylim(0,0.20)+
theme_bw() ->
r2.line
ggsave(r2.line, 
file = "r2.line.pdf",
w = 5, h = 2.5)

 #### --- plot nr99
ggplot() +
  geom_rect(data = final.windows.pos, 
  aes(xmin=start/1e6 , xmax =end/1e6, 
  ymin = 0, ymax = 0.025),
              fill = "gold", alpha = 0.6) +
geom_vline(xintercept = 2225744/1e6, color = "red") +
geom_vline(xintercept = 13154180/1e6, color = "red") +
geom_line(
data = win.out,
aes(
x=pos/1e6,
y=r2.99n/nloc
)) + 
ylab("proportion of r2 > 99") +
xlab("Chr 2L Mbp") +
#ylim(0,0.20)+
theme_bw() ->
r2.99n.line
ggsave(r2.99n.line, 
file = "r2.99n.line.pdf",
w = 5, h = 2.5) 
 
ggplot() +
  geom_rect(data = final.windows.pos, 
  aes(xmin=start/1e6 , xmax =end/1e6, 
  ymin = 0, ymax = 0.08),
              fill = "gold", alpha = 0.6) +
geom_vline(xintercept = 2225744/1e6, color = "red") +
geom_vline(xintercept = 13154180/1e6, color = "red") +
geom_line(
data = win.out,
aes(
x=pos/1e6,
y=r2.70n/nloc
)) + 
ylab("proportion of r2 > 70") +
xlab("Chr 2L Mbp") +
#ylim(0,0.20)+
theme_bw() ->
r2.70n.line
ggsave(r2.70n.line, 
file = "r2.70n.line.pdf",
w = 5, h = 2.5)

#####
   genofile <- seqOpen("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)
  ### common SNP.dt
    load("/netfiles/nunezlab/Drosophila_resources/Filtering_files/snp_dt_25percMissing.Rdata")
    setkey(snp.dt, id)

samps <- fread("/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/DEST_10Mar2021_POP_metadata.csv")


    getData <- function(chr="2L", start=14617051, end=14617051, samplesUse=samps$sampleId) {
      # chr="2L"; start=14617051; end=14617051; samplesUse=samps$sampleId

      ### filter to target
        snp.tmp <- data.table(chr=chr, pos=start:end)
        setkey(snp.tmp, chr, pos)
        setkey(snp.dt, chr, pos)
        seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id, sample.id=samplesUse, verbose=T)
    ### get annotations
      message("Annotations")
      tmp <- seqGetData(genofile, "annotation/info/ANN")
      len1 <- tmp$length
      len2 <- tmp$data

      snp.dt1 <- data.table(len=rep(len1, times=len1),
                            ann=len2,
                            id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))

    # Extract data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
      snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                            list(variant.id=id)]

      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
      snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]
      
      ### get frequencies
        message("Allele Freqs")

        ad <- seqGetData(genofile, "annotation/format/AD")
        dp <- seqGetData(genofile, "annotation/format/DP")

        if(class(dp)[1]!="SeqVarDataList") {
          dp.list <- list()
          dp.list$data <- dp
          dp <- dp.list
        }

        af <- data.table(ad=expand.grid(ad$data)[,1],
                         dp=expand.grid(dp$data)[,1],
                         sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                         variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

      ### tack them together
        message("merge")
        #afi <- merge(af, snp.dt1.an, by="variant.id")
        afi <- merge(af, snp.dt, by.x="variant.id", by.y="id")

        afi[,af:=ad/dp]

      ### calculate effective read-depth
        afis <- merge(afi, samps, by="sampleId")

        afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
        afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
        afis[,af_nEff:=round(af*nEff)/nEff]

      ### return
        #afis[,-c("n"), with=F]
        cbind(afis, snp.dt1.an)

    }

getData(chr = "2L", start = 5155959, end = 5155959  )

