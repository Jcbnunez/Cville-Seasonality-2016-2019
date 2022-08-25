# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  library(tidyr)

### load core20 GLM object
  load(file="/project/berglandlab/alan/core20glm.Rdata")

### load thermal GLM object
  job <- 0
  glm.fn <- paste("/project/berglandlab/thermal_glm_dest/processedGLM/glm.out.VA_ch_", job, ".Rdata", sep="")
  load(glm.fn)

### merge
  m1 <- merge(glm.out[mod=="aveTemp+year_factor"], core20.glm[mod=="season+locality_factor"][set=="machado"], by.x="id", by.y="variant.id")
  m2 <- merge(glm.out[mod=="aveTemp+year_factor"], core20.glm[mod=="season+locality_factor"][set=="no_va"], by.x="id", by.y="variant.id")

### clean up
  rm(glm.out, core20.glm)

### merge with haplotags
  load("/project/berglandlab/jcbnunez/Shared_w_Alan/haplo_tags_SNPids.Rdata")
  haplo_tags_SNPids_and_inv <- as.data.table(haplo_tags_SNPids_and_inv)
  haplo_tags_SNPids_and_inv[,pos:=as.numeric(pos)]
  setkey(m1, chr, pos)
  setkey(m2, chr, pos)

  setkey(haplo_tags_SNPids_and_inv, chr, pos)

  m1 <- merge(m1, haplo_tags_SNPids_and_inv, all.x=T)
  m2 <- merge(m2, haplo_tags_SNPids_and_inv, all.x=T)

### rank normalization
  m1 <- m1[chr!="X"]
  m2 <- m2[chr!="X"]

  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]


  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]


### get betas
  m1.beta <- m1[,list(thermal.beta=tstrsplit(b.x, ";")%>%last%>%as.numeric, seas.beta=tstrsplit(b.y, ";")%>%last%>%as.numeric), list(id)]
  m2.beta <- m2[,list(thermal.beta=tstrsplit(b.x, ";")%>%last%>%as.numeric, seas.beta=tstrsplit(b.y, ";")%>%last%>%as.numeric), list(id)]

  m1 <- merge(m1, m1.beta, by="id")
  m2 <- merge(m2, m2.beta, by="id")


### enrichment test
  m1[,inv:=invName!="none"]
  m2[,inv:=invName!="none"]

### save
  save(m1, m2, file="/scratch/aob2x/core20_dest_aveTemp.Rdata")

### load
  load(file="/scratch/aob2x/core20_dest_aveTemp.Rdata")


### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  setkey(m2, "chr")
  wins <- foreach(chr.i=c("2L", "2R", "3L", "3R"), .combine="rbind", .errorhandling="remove")%dopar%{
      tmp <- m2[J(chr.i)]
      data.table(chr=chr.i,
                  start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                  end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
  }
  wins[,i:=1:dim(wins)[1]]
  dim(wins)


### iterate through
  setkey(m2, chr, pos)
  thr.i <- 0.1

  win.out <- foreach(win.i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
      # win.i <- 50
      message(paste(win.i, dim(wins)[1], sep=" / "))
      win.tmp <- m2[J(data.table(chr=wins[win.i]$chr, pos=wins[win.i]$start:wins[win.i]$end, key="chr,pos")), nomatch=0]

     ### DEST set, enrichment
       tab <- table(win.tmp$thermal.rank <thr.i,
                    win.tmp$season.rank  <thr.i)

       fet <- fisher.test(tab)

     ### DEST set, sign test
       st.T <- sum(sign( win.tmp[thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) ==
                    sign(win.tmp[thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
       st.F <- sum(sign( win.tmp[thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) !=
                    sign(win.tmp[thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

       bt <- binom.test(st.T, st.T+st.F, .5)


       tmp2 <- data.table(win.i=win.i,
                  mod="core20_DEST", set="aveTemp_yearFactor",
                  or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                  st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])

      ### return
        return(tmp2)
  }


  win.out <- rbindlist(win.out)
  win.out <- merge(win.out, wins, by.x="win.i", by.y="i", all=T)



  save(win.out, file="~/core20_enrichment_window.Rdata")

### scp
  scp aob2x@rivanna.hpc.virginia.edu:~/core20_enrichment_window.Rdata ~/.

### library
  library(data.table)
  library(ggplot2)
  library(patchwork)
  load("~/core20_enrichment_window.Rdata")


  win.ou

  en <-
  ggplot(data=win.out, aes(x=log2(or), y=-log10(p), color=chr)) + geom_point()

  st <-
  ggplot(data=win.out, aes(x=st.pr, y=-log10(p), color=chr)) + geom_point()

  en + st


  ggplot(data=win.out, aes(x=start/2+end/2, y=sign(st.pr-.5)*log10(p), color=or)) +
  geom_point() +
  facet_grid(~chr)








  thrs <- expand.grid(sapply(1:9, function(x) x*10^(-5:-1)))[,1]
  #thrs <- c(0.01, 0.05)

  o <- foreach(chr.i=unique(m1$chr), .combine="rbind", .errorhandling="remove")%do%{
    foreach(inv.i=c(T,F), .combine="rbind", .errorhandling="remove")%do%{
      foreach(thr.i=thrs, .combine="rbind", .errorhandling="remove")%dopar%{
        # chr.i <- "2L"; inv.i <- T; thr.i<-0.05
          message(paste(chr.i, inv.i, thr.i, sep=" / "))

          ### Machado set, enrichment
            tab <- table(m1[chr==chr.i][inv==inv.i]$thermal.rank  <thr.i,
                         m1[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

            fet <- fisher.test(tab)

         ### Machado set, sign test
           st.T <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) ==
                        sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta.y))
           st.F <- sum(sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) !=
                        sign(m1[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta.y))

           bt <- binom.test(st.T, st.T+st.F, .5)


           tmp1 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i,
                      mod="core20_Machado", set="aveTemp_yearFactor",
                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])


         ### DEST set, enrichment
           tab <- table(m2[chr==chr.i][inv==inv.i]$thermal.rank  <thr.i,
                        m2[chr==chr.i][inv==inv.i]$season.rank <  thr.i)

           fet <- fisher.test(tab)

         ### DEST set, sign test
           st.T <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) ==
                        sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))
           st.F <- sum(sign( m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$thermal.beta) !=
                        sign(m2[chr==chr.i][inv==inv.i][thermal.rank<thr.i & season.rank<thr.i]$seas.beta))

           bt <- binom.test(st.T, st.T+st.F, .5)


           tmp2 <- data.table(chr=chr.i, inv=inv.i, thr=thr.i,
                      mod="core20_DEST", set="aveTemp_yearFactor",
                      or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2],
                      st.pr=bt$estimate, st.p=bt$p.value, st.lci=bt$conf.int[1], st.uci=bt$conf.int[2])

        ### return
          rbind(tmp1, tmp2)

      }
    }
  }

  o.dest_core20 <- o


### WINDOW - haplotag SNPs
  m1[season.rank<.05,list(pr=mean(sign(thermal.beta)==sign(seas.beta.y), na.rm=T), .N), list(win)]
  m2[season.rank<.05,list(pr=mean(sign(thermal.beta)==sign(seas.beta), na.rm=T), .N), list(win)]

  m2[,list(TT=sum(season.rank<.05 & thermal.rank<.05),
          TF=sum(season.rank<.05 & thermal.rank>.05),
           FT=sum(season.rank>.05 & thermal.rank<.05),
          FF=sum(season.rank>.05 & thermal.rank>.05)), list(win)]


  save(o.dest_core20, file="~/dest_core20.enrichment_signTest.Rdata")















  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m1[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]


  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,thermal.rank := rank(p.lrt.x)/length(p.lrt.x)]
  m2[!is.na(p.lrt.x) & !is.na(p.lrt.y) ,season.rank := rank(p.lrt.y)/length(p.lrt.y)]
