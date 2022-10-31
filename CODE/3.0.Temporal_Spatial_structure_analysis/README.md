# 3.0.Temporal_Spatial_structure_analysis

### 1.Import_GDStoR.r
This code import the `gds` object into R. Generates allele frequency files, filters for effective coverage and outputs R objects for each chromosome.

### 1.1.SBATCH_import_GDS_R.sh
This is the batch launcher of 1.Import_GDStoR.r

### 1.3.Rdat_objects_ECfiltered
This is a folder with R objects from 1.Import_GDStoR.r. These objects are filtered.

### 2.make_AllDEst_Robject.R
Do some additional filtering and prepare a single R object for the multivariate analysis

### 3.Make_PCA_time_space_DEST.R
Conducts PCA on the data object

### 3.1.SBATCH_make_DEST_PCA.sh
launcher of 3.Make_PCA_time_space_DEST.R

### 4.Temporal_Spatial_Analysis_ofPCA.R
Run analyses on the output of the PCA

### 5.Calculate_Fst_DEST_global.R
Calculated FST on the gds obejct (overwintering FSTs)

### 6.Calculate_Fst_CVILLE.window.R
An implementation of FST in windows for Charlottesville

### 6.1.launch.winFSTCville.sh 
Launcher of 6.Calculate_Fst_CVILLE.window.R

## Folders 
Contain metadata and some outputs for this analysis
* dest_filtering
* Figures
* guide_files
* methods_figures
* tables