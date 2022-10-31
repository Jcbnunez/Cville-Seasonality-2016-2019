### Organize files for mapping
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(foreach)

setwd("/scratch/yey2sn/Overwintering_ms/msp300.case/fasta.files.fld")

sech.met <- vroom("/scratch/yey2sn/Overwintering_ms/msp300.case/Sechelia.metadata.txt")
yak.met <- vroom("/scratch/yey2sn/Overwintering_ms/msp300.case/yakuba.metadata.txt")

gen.metadat = function(dat, sp){
  metadat = foreach(i=1:dim(dat)[1],
                    .combine = "rbind")%do%{
                      message(i)
                      
                      grep.files <- system(paste("ls | grep .sra"), intern = T )
                      readsof.int <- grep.files[grep(dat$SRR[i], grep.files)]
                      
                      grep.f <- system(paste("ls | grep ", dat$SRR[i]), intern = T )
                      
                      n.read <- length(grep.f)
                      
                      if(n.read == 1){
                        message("single end")
                        data.frame(
                          sampId = dat$samp[i],
                          config = n.read,
                          sp = sp,
                          SRR = dat$SRR[i],
                          readF = grep.f,
                          readR = NA
                        )
                      } else
                        if(n.read == 2){
                          message("pair end")
                          data.frame(
                            sampId = dat$samp[i],
                            config = n.read,
                            sp = sp,
                            SRR = dat$SRR[i],
                            readF = grep.f[1],
                            readR = grep.f[2]
                          )
                        }
                      
                    }
  
  write.table(metadat, 
              file = paste(sp,".metadat.mapping.txt", sep = ""), 
              append = FALSE, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
}

gen.metadat(sech.met, "sech")
gen.metadat(yak.met, "yak")
