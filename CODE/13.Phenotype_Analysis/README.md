# 13.Phenotype_Analysis

## GWAS

### 1.GWAS_DIRS.R
This script creates directories to 1. hold gwas files, 2: creates a reference table read by gwas scripts, 3: creates permuted phenotype data files for use in permutation gwas creation. 

#### Files needed:
* `updatedphenos`: A list of the phenotypes used in the analysis with their respective line average for the DGRP. This can be found as a supplementary table in the main paper. **This is Table S11** when downloaded as saved as "updatedphenos.txt"**
* This code creates: sets of permutations for selected phenos

### 2.GWAS_WithGRMS.R
Creating a GWAS file and do the GWAS. This code uses a GRM.

#### Files needed:
* `dgrp.gds`: The DGRP data in GDS format
* `new.gwas.ref.table.txt`: created in 1.
* The inversion status of the DGRP lines, found the DGRP website. 
* The wolbachia status file, also found the DGRP website. 

### 3.GWAS_WithoutGRMS.R
Creating a GWAS file and do the GWAS. This code does not uses a GRM. Instead uses an identity matrix. 

## Permutation code of GWAS
These permutations are not featured in the paper. but we provide the code for those interested.

### 4.GWAS_WithGRM_perms .R
Identical to 2 but uses permuted data

### 5.GWAS_WithoutGRM_perms.R
Identical to 3 but uses permuted data

### 6.GMMATdatabin.R
Combines data 2,3,4, and 5.

### 7.PhenoStats.R
Compares GWAS with and w/o GRM.

## Inversion - Phenotype association
These bits of code run explicit association analyses between the phenotype data and the inversion status of the lines. 

### 8.Inversion.modeling.comp.R
Make linear models looking at the relationship between line status and phenotype value.

### 9.crescent.analysis.R
This does an enrichment between the GWAS analysis (see code 3; no GRM) and the seasonal GLM analysis. looking for phenotypes with SNPs that are enriched in both analyses. 

```
#!/usr/bin/env bash
#
###goal, creat a slurm script that permutes gmmat for different phenotypes
#SBATCH -J pheno_glm_enrichment # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH --mail-type=end
#SBATCH --mail-user=bal7cg@virginia.edu
#SBATCH -o /scratch/bal7cg/score_error/surv.gmmat.%A_%a.err # Standard error
#SBATCH -e /scratch/bal7cg/score_output/surv.gmmat.%A_%a.out # Standard output
#SBATCH -p largemem
#SBATCH -A berglandlab
#SBATCH --array=1-100


module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0
module load intel/18.0 intelmpi/18.0

#define variables

jobid=${SLURM_ARRAY_TASK_ID}

Rscript \
--vanilla \
11.crescent.analysis.R \
${jobid} \
```

### 10.crescentplot.filtering.R
This collects, filters (i.e., permutation analysis), and summarizes the output of 9.

### 11.sliding.window.jobfile.R
Makes a job file to be use for 13. This job is used in an array job using a SLURM manager:
```
#!/usr/bin/env bash
#
###goal, creat a slurm script that gathers data from permutation files
#SBATCH -J pheno_gather # A single job name for the array
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 #<= this may depend on your resources
#SBATCH --mem 120G #<= this may depend on your resources
#SBATCH --mail-type=end
#SBATCH --mail-user=something@virginia.edu
#SBATCH -o ./score_error/surv.gmmat.%A_%a.err # Standard error
#SBATCH -e ./score_output/surv.gmmat.%A_%a.out # Standard output
#SBATCH -p standard
#SBATCH -A berglandlab
#SBATCH --array=2-127%10

module load goolf/7.1.0_3.1.4
module load gdal proj R/4.0.0
module load intel/18.0 intelmpi/18.0

#define variables

jobid=${SLURM_ARRAY_TASK_ID}
wd="/project/berglandlab/Yang_Adam/" #working directory

Rscript \
--vanilla \
${wd}Pheno_gather_Sep_14.R \
${jobid} \
```

### 12.slidingwindowanalysis.R
Similar to 9, but we implemented across sliding windows.

### 13.slidingwindowdata.filter.R
Collect and summarize output of 12. Filter by Bonferroni and permutations

### 14.pca.phenoplot.R
Filter down phenotyes for 13. Conduct PCA on these data and downsetram analyses.

## Other Folders

### DEPRECATED

