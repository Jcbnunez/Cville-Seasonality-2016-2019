#### Create square matrix from DF in R

library(data.table)
library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
file_in=args[1]
window_in=args[2]
window_file=args[3]

#file_in="./output/win.103.ld"
#window_in="win.103"
#window_file= "./windows/win.103.txt"

message(file_in) 
message(window_in) 
message(window_file) 

##read in file
metadat_tmp = fread(window_file, header = F)
metadat_tmp %<>% separate(V1, into = c("chr", "pos", "type"), sep = "_" )

in_tmp = fread(file_in)
in_tmp %>%
  filter(BP_A %in% metadat_tmp$pos  ) %>%
  filter(BP_B %in% metadat_tmp$pos) ->
  in_tmp_metadat

in_tmp_metadat$BP_A %>% unique %>% length
in_tmp_metadat$BP_B %>% unique %>% length

### Make square matrix
in_tmp_metadat %>%
  dcast(SNP_A~SNP_B, value.var = "R2") ->
  square_matrix_ld

### save file
final_file_destin="./square_matrices_final"

save(square_matrix_ld, file= paste(final_file_destin, "/" , window_in, ".squareMatrix.LD.Rdata", sep = ""))

