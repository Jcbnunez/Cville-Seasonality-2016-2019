### Make new figure 3
### 

### libraries
library(tidyverse)
library(data.table)
library(tidyr)
library(viridis)
library(patchwork)
library(magrittr)
library(RColorBrewer)
library(rstatix)

##### input
output_results_window <- "/scratch/yey2sn/Overwintering_ms/4.1.NewTempModels/all.mod.out.Rdata"

inversion_map <- "/project/berglandlab/Dmel_genomic_resources/Inversions/InversionsMap_hglft_v6_inv_startStop.txt"

####
####
load(output_results_window)


### summarize and plot
all.mod.out %>%
  mutate(model_name = paste(label,start,end, sep ="_" )) %>%
  group_by(model_name,chr,invName,  perm_type ) %>%
  summarize(med_rnvp = qlogis(median(rnp.binom.p)),
            rnvp90 = qlogis(quantile(rnp.binom.p,0.9)),
            rnvp10 = qlogis(quantile(rnp.binom.p,0.1)),
  )

