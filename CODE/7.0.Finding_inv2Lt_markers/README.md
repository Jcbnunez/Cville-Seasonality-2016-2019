# 7.0.Finding_inv2Lt_markers

This folder contains code to generate our inversion, In(2L)t, markers. Here is code to also train a SVM to predict inversion status.


### 1.Process_data_from_vcf.R
Extract data from the DGRP vcf.

### 2.1.launch_make_pca.sh
Launcher of  2.make_first_pca.r

### 2.make_first_pca.r
Makes a PCA with the DGRP data. The PC signal in 2L is driven by inversion status.

### 3.run_pca_correlations.r
Run correlation analysis of individual SNP contributions to PCA. Preliminary markers are created here. 

### 3.1.launch_run_pca_correlations.sh
Launcher of 3.run_pca_correlations.r

### 4.mine_correlations_to_PCA.r
Examine the output of 3.run_pca_correlations.r

### 5.calculate_ld_vcftools.sh
Calculate LD among markers using VCF. Another strategy using plink is shown below. In the DGRP

### 6.1.plot_ld_from_DGRP.R
Visualize LD In the DGRP

### 6.calculate_ld_plink.sh
Another strategy to calculate LD, using plink. In the DGRP
```
#Example -- rank of LD 0.99-0.90
sbatch --array=1-750 \
6.calculate_ld_plink.sh \
Rank90_99
```
### 6.launch_mine_ld_inR.sh
Launcher of  6.mine_ld_data.r

### 6.mine_ld_data.r
Parse and summarize the plink output

### 7.liftOver.inversionMarkers.sh -- main version
Lift Over across genomes using the stand-alone program LiftOver. In the DGRP

### 7.prepare_snpids_ldcalcc.r -- alt. version
Lift Over across genomes using R, a version with the bash program is shown above. In the DGRP

### 7.make.2lt.bed.r
Make in2Lt bed file

### 7.translate_dm3_2ltsnos_to_dm6.r
Noodle around with lifover data and consolidate into one object

### 8.calculate_ld_plink_inCM.sh
Use plink to calculate LD in the whole genome data from wild flies collected in Virginia.
```
#0.99-0.90
sbatch --array=1-77 \
8.calculate_ld_plink_inCM.sh \
Rank99
```
### 8.mine_ld_dat_CM.r
Same as in the DGRP. Mine and summarize VA LD data.

### 8.sanity_check_dgrp_CM.r
A sanity check ...

### 9.create_vcf_for_dapc_pred.sh
Create a VCF file for an exploratory DAPC analysis

### 10.train_predictive_model.r
Train the SVM model. Save it to an R object for later deployment. 
