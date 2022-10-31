# 11.0.FST_inversion_pooled

This folder contains code for $F_{ST}$ analyses done using [`poolfstat`](https://cran.r-project.org/web/packages/poolfstat/index.html). These analyses look only at markers in the inversion In(2L)t.  There are two routes to this analysis and both are used in the paper.

## Pop-wise
The pop-wise pipeline calculates mean $F_{ST}$ within and between population sets: Charlottesville, EU-E, EU-W. Code is organized by population set inside the folders

## SNP-wise
The SNP-wise pipeline calculates mean $F_{ST}$ within and between population sets constrained for particular loci in the genome. This is a computationally demanding pipeline. Code is named based on its annotation and should be run in order. 