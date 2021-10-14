### consolidate the bk ld output  file
### in R

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

library(data.table)
library(tidyverse)
library(magrittr)

#read in all the files
plink_out <- system("ls bk_ld_output", intern = T)
root <- "/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/bk_ld_output/"

out_list = list()
for(i in 1:length(plink_out)){

  tmp <- fread( paste(root,plink_out[i], sep = "" ) )
  
  if( dim(tmp)[2] == 6 ) {
    print("ok")
    
    tmp %<>%
      dcast(V2+V1~V6, value.var = "V4") %>% 
      separate(V2, into = c("focal_snp", "affected_snp"), sep ="vs")
    
    out_list[[i]] = tmp
    
  } else {print("does not comply")}
  
}

bk_ld_frequencies = do.call(rbind, out_list)
  
save(bk_ld_frequencies, file = "bk_ld_frequencies.Rdata")  


