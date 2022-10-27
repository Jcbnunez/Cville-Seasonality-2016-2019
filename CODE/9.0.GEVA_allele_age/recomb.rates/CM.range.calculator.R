#### CM calcualtor
#### 

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)


## based on: https://www.sciencedirect.com/science/article/pii/S0378111910001769
## 
CM.calc.2L <- function(bp){
  
  x = bp/1e6
  x = ifelse(x > 18.1, 18.1, x)
  
  cm <- -0.01*(x^3)+0.20*(x^2)+2.59*(x)-1.56
  
  cm = ifelse(cm < 0, 0, cm)

  return(cm)
  
}

cm.plot = foreach(i=seq(from = 0, to = 25e6, by = 1e5), 
                  .combine = "rbind")%do%{
  data.frame(pos=i, cm=CM.calc.2L(i))
}

ggplot(data = cm.plot, aes(x=pos, y =cm )) + geom_line()

#### Add recombination map info to map file
map.file.come <- fread("/Users/jcbnunez/Documents/GitHub/Cville-Seasonality-2016-2019/9.GEVA_allele_age/recomb.rates/chr2L.comeron.RcRate.txt", header = T)

###estimate CM
map.file.come %>%
  group_by(start,end,chr) %>%
  mutate(CM.w = (CM.calc.2L(start) + CM.calc.2L(end))/2) %>%
  mutate(med.pos =median(start:end) ) %>%
  dplyr::select(
    "Chromosome"=chr,
    "Position(bp)"=med.pos,
    "Rate(cM/Mb)"=Comeron_cm_mb,
    "Map(cM)"=CM.w
  ) %>%
  .[,-c(1:2)] ->
  map.file.2L.GEVA

write.table(map.file.2L.GEVA, file = "map.file.2L.GEVA.map", 
            append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

