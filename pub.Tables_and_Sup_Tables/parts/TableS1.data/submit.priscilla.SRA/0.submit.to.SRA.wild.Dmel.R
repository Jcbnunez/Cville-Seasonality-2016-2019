### submit prsicilla SRA.
### 

library(tidyverse)
library(reshape2)
library(magrittr)
library(data.table)
library(foreach)

metadat <- fread("/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/samps.used.txt")
metadat %<>%
  separate(sampleId, remove = F, 
           into = c("pop", "well", "etc"),
           sep = "_"
           )
  
files.loc <- "/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data"

files.to.submit = foreach(i = 1:dim(metadat)[1],
                          .combine = "rbind")%do%{
  
  tmp.files <- list.files(paste(files.loc, paste(metadat$pop[i], metadat$well[i], sep = "_"), sep = "/"  ))
  
  tmp.files.reads <- tmp.files[grep("fq.gz",tmp.files)]
  
  forward.read = tmp.files.reads[grep("_1.", tmp.files.reads)]
  reverse.read = tmp.files.reads[grep("_2.", tmp.files.reads)]
  
  data.frame(pop = metadat$pop[i],
             well = metadat$well[i],
             read1 = forward.read,
             read2 = reverse.read,
             address1=paste(paste(files.loc, paste(metadat$pop[i], metadat$well[i], sep = "_"), sep = "/"  ), forward.read, sep = "/"),
             address2=paste(paste(files.loc, paste(metadat$pop[i], metadat$well[i], sep = "_"), sep = "/"  ), reverse.read, sep = "/")
  )
    
}


full_join(metadat, files.to.submit) ->
  master.sra.file.OW

write.table(master.sra.file.OW, 
            file = "master.sra.file.OW.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
