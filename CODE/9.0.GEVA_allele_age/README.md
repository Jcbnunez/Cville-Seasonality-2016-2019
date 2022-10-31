# 9.0.GEVA_allele_age
Lear more about GEVA here: https://github.com/pkalbers/geva 

### 1.prepare.analysis.VCFs.sh
A script that parses VCF files in preparation for the GEVA TMRCA analysis

### 2.run_geva_dmel.sh
Master script that runs the program GEVA
```
sbatch \
--array=1-$( cat array_job_guidefile.all.2L.only.txt | wc -l) \
2.run_geva_dmel.sh \
array_job_guidefile.all.2L.only.txt \
0.01 \
Dmel_ALL \
CM_pops.2L.all.phased.recode.vcf.gz

sbatch \
--array=1-$( cat array_job_guidefile.all.2L.only.txt | wc -l) \
2.run_geva_dmel.sh \
array_job_guidefile.all.2L.only.txt \
0.01 \
Dmel_het \
CM_pops.2L.het.phased.recode.vcf.gz

sbatch \
--array=1-$( cat array_job_guidefile.all.2L.only.txt | wc -l) \
2.run_geva_dmel.sh \
array_job_guidefile.all.2L.only.txt \
0.01 \
Dmel_inv \
CM_pops.2L.inv.phased.recode.vcf.gz

sbatch \
--array=1-$( cat array_job_guidefile.all.2L.only.txt | wc -l) \
2.run_geva_dmel.sh \
array_job_guidefile.all.2L.only.txt \
0.01 \
Dmel_std \
CM_pops.2L.std.phased.recode.vcf.gz

```

### 3.explore.geva.results.R
Explore results in R

## Extra
* recomb.rates: This file contains recombination rate information from melanogaster. See https://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl
* array_job_guides
* chromosome_name_file
* extra
* old_deprecated
* sample_name_guides
