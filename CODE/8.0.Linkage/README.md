# 8. Linkage

### 1.Prepare_SNP_data.R
Extract the SNP ids for the inversion that are outlier SNPs in the GLM -- In2Lt
 
### 2.make_vcf_for_plink.sh
Prepare a VCF, using VA data, that has the desired markers for input into Plink.

### 3.run_plink.sh
Calculate LD in plink

### 4.1.collect_and_bind_ld_dopar.R
run an R script that collects and summarizes the LD data

### 4.2.run_collect_ld_dopar.sh
Launcher of 4.1.collect_and_bind_ld_dopar.R

### 5.parse_LD_data.R
Parse and noodle around with the LD data output

### 6.pull_annotation_GLM_LD_snps.R
Annotate the identity of each mutation of interest.

### 7.graph_ld_relationships.R
Make some preliminary graphs of the LD relationships.

### 8.matched_controls_survey.R
This code identifies a series of matched controls relative to GLM markers and estimates LD for those. This analysis is deprecated because it is not part of the paper, but the code is here for those interested in doing this, need be. These sub-analysis have the following code bits:

* 9.make_vcf_for_plink_matchedControls.sh
* 10.run_plink_matched_controls.sh
* 11.collect_and_bind_ld_dopar_CONTROL.R
* 11.run_collect_ld_dopar_CONTROL.sh
* 12.parse_LD_data.control.R


### 13.1.make_vcf_for_plink_intrahap.sh
This is another bit of aux code. This partitions the VCF into all inverted homozygous samples and all standard homozygous samples. For LD purposes. 

### 13.2.run_plink_intra_haplotype.INV.CM.sh
run plink LD calcl for the homozygous samples sets. The follow-up steps can be found here:

* 14.1.collect_and_bind_ld_INVhaps.R
* 15.graph_ld_relationships.INV.R

### 16.PLOT_LD_FIGURE.final.R
This bit of code collects the appropriate data generated above and produces a figure.

## Folders

* snp_guide_files
* individual_guide_files
* DEPRECATED
* Figures
* outputs
