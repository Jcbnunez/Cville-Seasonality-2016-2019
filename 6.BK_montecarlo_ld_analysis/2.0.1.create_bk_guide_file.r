### make guide files for the iterators
### in R

#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

library(data.table)
library(tidyverse)
library(magrittr)

# files to read:
glm_out <- "./final_in2Lt_markers.txt"
inv_wlt_out <- "./temperature_snps_ids_top100.txt"
#

# files to iterate:
iterate_over_inv="/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/final_in2Lt_markers.txt"
iterate_over_temp="/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/temperature_snps_ids_top100.txt"

###
glm_dt <- fread(glm_out, header = FALSE)
glm_dt %<>%
	mutate(type = "temperature",
			iterate_file = iterate_over_inv)


inv_ft <- fread(inv_wlt_out, header = FALSE)
inv_ft %<>%
	mutate(type = "inversion",
			iterate_file = iterate_over_temp)


rbind(glm_dt, inv_ft) -> guide_files_bk_ld

write.table(guide_files_bk_ld,
            file = "./guide_files_bk_ld.txt",
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)