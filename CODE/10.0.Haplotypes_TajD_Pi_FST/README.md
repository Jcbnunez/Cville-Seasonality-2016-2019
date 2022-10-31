# 10.0.Haplotypes_TajD_Pi_FST
This folder contains code for a battery of tests and visualizations related to haplotype strcuture and site frequency spectrum analyses. 

### 1.Haplotype_analyses.R
Visualizes haplotypes in the data, and generates files with markers of interest.

### 3.handle_dna_seq.R
Uses an MSA reader in R to import fasta sequences 

### 5.multidimentional_analysis.R
Runs PCA on the inversion data.

* 5.2.launch_dimdesc.sh: Launcher file
* 5.1.run.dim.desc.r: run the correlation analysis of the PC projects

### 6.explore_pca_dimdesc.r
Explores and summarizes the output of the PCA analysis

### 7.1.prepare_vcf_for_plotting.sh
Prepare the joint VCF files for plotting 

### 7.2.process_haplotypes.R
Uses R to process the VCF files. Extract inversion and GLM markers in haplotypes

### 7.3.plot_haplotypes.R
Plots data

### 8.2.1.makeBEDfile.R
Makes bed file with inversion and GLM data

### 8.2.2.make.markers.for.fst.test.R
Extracts markers of interest for a downstream FST test

### 8.2.run_hapCount_CM_subdivided.sh
Run the VCFTools function `hapcount`

### 8.2.run_pi_D_CM_subdivided.sh
Run $\pi$,  Tajimas D stats on the data

### 8.2.run_pi_D_CM_for5.1_subdivided.sh
Run $\pi$,  Tajimas D stats on the data, for window w5.2 (w = 5.1-5.2)

### 9.1.collect_vcftools_analyses.R
A script which collects the output of various VCFtools analyses and summarizes them

### 8.3.run_fst.sh
Run $F_{ST}$ on the data using VCFtools

### 9.4.Fst_DGRPvsCM.plot.R
plot $F_{ST}$ in the Virginia data and the DGRP data

### 10.create_makers_for_window_haps.R
This script creates a series of marker files in preparation for the big haplotype plot shown in the paper

### 11.extract_vcf_for_haploplotting.sh
This script extracts the SNP markers for the  big haplotype plot shown in the paper

### 12.1.plot_haplotype.homozyg.PCA.windows.R
Plot the big haplotype plot shown in the paper

### 14.Code.SpecialFor.DPGP_FST
This is code for calculating $F_{ST}$ in the DPGP3.

## Plotting scrtips

* 9.2.Plot_FST_CM_INVvsSTD_analyses.R
*  9.2.Plot_vcftools_analyses.For5.1.R
* 9.2.Plot_vcftools_analyses.R
* 9.3.plot_haplo.counts.R

## Supplementary

### 13.select_dgrp_lines_for_adam.R
Used for deficiency mapping downstream

## Extra
* guides: Contains important guide files for analyses above
* data: some small data files used here
* Deficiency_mapping_Adam_lines
* DEPRECATED
* Figures