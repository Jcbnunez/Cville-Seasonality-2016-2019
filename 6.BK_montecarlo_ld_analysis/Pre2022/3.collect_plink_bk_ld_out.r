### consolidate the bk ld output  file
### in R

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

library(data.table)
library(tidyverse)
library(magrittr)

# load metadata of comparison
met_f <- "/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/guide_files_bk_ld.txt"
meta_dat <- fread(met_f, header = F)
meta_dat %<>%
  mutate(comparison = paste(V1, V3, sep = "|"))
names(meta_dat) = c("snp_1","snp_1_type", "snp_2","snp_2_type", "test_type" , "comparison")

#read in all the files
plink_out <- system("ls bk_ld_output", intern = T)
root <- "/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/bk_ld_output/"

out_list = list()
for(i in 1:length(plink_out)){
  
  tmp <- fread( paste(root,plink_out[i], sep = "" ), header = F )
  
  if( dim(tmp)[2] == 5 ) {
    print("ok")
    
    names(tmp) = c("comparison", "gametic_phase", "Freq", "Freq_LE", "status")
    
    tmp %<>%
      dcast(comparison~status, value.var = "Freq") %>%
      left_join( meta_dat)
    
    out_list[[i]] = tmp
    
  } else {print("does not comply")}
  
}

bk_ld_frequencies = do.call(rbind, out_list)
  
save(bk_ld_frequencies, file = "bk_ld_frequencies.Rdata")  

