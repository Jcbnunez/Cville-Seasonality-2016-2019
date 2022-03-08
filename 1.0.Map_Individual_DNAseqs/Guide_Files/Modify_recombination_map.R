library(tidyverse)

hglft_v6_recombination_map <- read.delim("~/Dropbox/2021_Cville_DEST/Alyssa's Single Ind/Guide_Files/hglft_v6_recombination_map.bed", header=FALSE)

hglft_v6_recombination_map %>%
  group_by(V1) %>%
  mutate(cM= cumsum(V4),
         "Pos(bp)"= round( (V2+V3)/2,0) ) ->
  Chr_recomb



names(Chr_recomb) = c("Chromosome","Start","End","Rate(cM/Mb)","Map(cM)","Position(bp)")
Chr_recomb$Chromosome = gsub("chr","",Chr_recomb$Chromosome)


write.table(file = "~/Dropbox/2021_Cville_DEST/Alyssa's Single Ind/Guide_Files/hglft_v6_GEVA_map.map",
            Chr_recomb[,c("Chromosome","Position(bp)","Rate(cM/Mb)","Map(cM)")],
           sep = "\t",
           row.names = F,
           col.names = T,
           quote = F
           )