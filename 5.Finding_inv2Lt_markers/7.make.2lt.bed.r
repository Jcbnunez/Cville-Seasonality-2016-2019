#Make dm3 inversion markers

#ijob -A jcbnunez -c 10 --mem=40G  --partition=standard
#module load intel/18.0 intelmpi/18.0
#module load goolf/7.1.0_3.1.4
#module load gdal proj R/4.0.0
#R

library(data.table)
library(tidyverse)
library(magrittr)

markers <- fread("inv2L_informative_markers_Dm3.txt")

markers_bed_pos <-
	data.frame(pos_format = 
			   paste(
			   paste("chr", markers$chr, sep = "" ),
			   ":",
			   markers$pos,
			   "-",
			   markers$pos,
			   sep = "")
	) 
	
write.table(markers_bed_pos, 
            file = "inv_2l_dm3_markers_asPos.txt", 
            append = F, 
            quote = F, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = F,
            col.names = F)

#as . bed
#

markers_bed_bed <-
  data.frame(pos_format = 
                 paste("chr", markers$chr, sep = "" ),
                 markers$pos-1,
                 markers$pos
  ) 

write.table(markers_bed_bed[order(markers_bed_bed[2]),], 
            file = "inv_2l_dm3_markers_asBed.bed", 
            append = F, 
            quote = F, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = F,
            col.names = F)
