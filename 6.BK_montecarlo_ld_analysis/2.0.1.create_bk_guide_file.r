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
# First lod in the inversion markers
### Note in the github this file can be found at: Cville-Seasonality-2016-2019/6.BK_montecarlo_ld_analysis/Files/
Inversion_markers <- "./final_in2Lt_markers.txt"
inv_dt <- fread(Inversion_markers, header = FALSE)

# Second, load in the glm markers
### Note in the github this file can be found at: Cville-Seasonality-2016-2019/6.BK_montecarlo_ld_analysis/Files/
glm_outliers <- "/scratch/yey2sn/Overwintering_ms/6.BK_test_montecarlo_sim/Cville_glm_p01_hits_annot.txt"
glm_dt <- fread(glm_outliers, header = FALSE)

#Create list of targets
c(inv_dt$V1, glm_dt$V43 ) -> loaded_markers

#validate markers in DEST
# In RIVANNA (UVA's supercomputer) this file can be found at: /project/berglandlab/jcbnunez/Shared_w_Alan
load("/scratch/yey2sn/Overwintering_ms/1.Make_Robjects_Analysis/Cville_2L.ECfiltered.Rdata")
SNPS_in_DEST <- paste( colnames(o), "SNP", sep = "_")

### validate SNPS
loaded_markers[loaded_markers %in% SNPS_in_DEST] -> valid_markers

data.frame(focal_snp = valid_markers) %>% 
  mutate(type= case_when(focal_snp %in% inv_dt$V1 ~ "inv_focus",
                         focal_snp %in% glm_dt$V43 ~ "tmp_focus",
                         )) -> validated_SNPs

## add extra metadata
validated_SNPs %<>%
  separate(focal_snp, 
           remove = F,
           into = c("chr", "pos", "feature"), sep = "_")
validated_SNPs$pos = as.numeric(validated_SNPs$pos)

validated_SNPs <- validated_SNPs[order(validated_SNPs$pos),]

#### thin by physical distance
#### This R script is housed in Rivanna. It is also contained in the Github of the project
source("/project/berglandlab/Dmel_genomic_resources/Software_and_Scripts/ThinLDinR_SNPtable.R")

validated_SNPs$chr = as.character(validated_SNPs$chr) 
validated_SNPs$pos = as.numeric(validated_SNPs$pos)

##

validated_SNPs[which(validated_SNPs$type == "inv_focus"),] -> inv_markers

validated_SNPs[which(validated_SNPs$type == "tmp_focus"),] -> tmp_markers

#thin inversion markers at 50K? -- optional, set to 1, 10, 100 ... if you do not want to thin
picksnps_inv <- pickSNPs(inv_markers,
                         dist=50000)

#thin glm .. or other markers at 50K? -- optional, set to 1, 10, 100 ... if you do not want to thin
picksnps_glm <- pickSNPs(tmp_markers,
                         dist=50000)

rbind(
inv_markers[picksnps_inv,],
tmp_markers[picksnps_glm,]) ->
  selected_markers_for_analysis

####
selected_markers_for_analysis %>%
  mutate(SNP_id = paste(focal_snp, type, sep = "|" )) %>%
  select(SNP_id) -> selected_markers_id

selected_markers_id$SNP_id = gsub("_focus", "", selected_markers_id$SNP_id)

expand.grid(selected_markers_id$SNP_id , selected_markers_id$SNP_id ) %>% 
  separate(Var1, 
           into = c("SNP_id_1", "type_1"),
           sep = "\\|") %>%
  separate(Var2, 
           into = c("SNP_id_2", "type_2"),
           sep = "\\|") %>% 
  mutate(comparison = paste(type_1,type_2, sep = ">") ) ->
  guide_files_bk_ld

##write.table(guide_files_bk_ld,
##            file = "./guide_files_bk_ld.txt",
##            append = FALSE, 
##            quote = FALSE, 
##            sep = "\t",
##            eol = "\n", 
##            na = "NA", 
##            dec = ".", 
##            row.names = FALSE,
##            col.names = FALSE)
##
### Make unique list
guide_files_bk_ld %>%
  .$SNP_id_1 %>%
  unique() ->
  unique_snp_list

write.table(unique_snp_list,
            file = "./unique_snp_guide_bk_ld.txt",
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            eol = "\n", 
            na = "NA", 
            dec = ".", 
            row.names = FALSE,
            col.names = FALSE)
