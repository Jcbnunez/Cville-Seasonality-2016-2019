#module load gcc/7.1.0
#module load openmpi/3.1.4
#module load R/4.1.1
#module load gdal
#module load proj

### in R
library(adegenet)

DPGP3 = fasta2genlight("./DPGP3.cat.fasta", 
quiet = FALSE, chunkSize = 10000, saveNbAlleles = FALSE,
               parallel = TRUE, n.cores = 2)
