# Repo for: "A cosmopolitan inversion drives seasonal adaptation in overwintering _Drosophila_"

This is the git repo for the paper "A cosmopolitan inversion drives seasonal adaptation in overwintering Drosophila". This repo contains all the code and data needed to replicate our analyses.

# Citation
Please cite our paper as TBD

## Analysis scheme flowchart
In our paper we conduct a series of analyses to show that the cosmopolitan inversion In(2L)t is a hot-spot of seasonal adaptation in temperate Drosophila. Our analyses use four general data types: 
1. First, pooled-seq data of seasonal flies collections. These data are a combination of new data generated for this paper as well as data from [DEST](https://dest.bio/). 
2. Second, individual whole genome data  from flies in Virginia. These data are also combines with other whole genome datasets of Drosophila, including the [DGRP2](http://dgrp2.gnets.ncsu.edu/), [DPGP3](https://www.johnpool.net/genomes.html), and others (see below). 
3. Third, we used published phenotype and GWAS data done on the [DGRP2](http://dgrp2.gnets.ncsu.edu/).
4. Complete gene sequences for the gene [_Msp300_](https://flybase.org/reports/FBgn0261836) in _Drosophila melanogaster_, _D. simulans_, _D. yakuba_, _D. sechelia_ and _D. mauritana_ 

Given the diversity of these data types we show the analyses workflows below and the corresponding script sets needed to reproduce the analyses. In these charts, squares represent data, diamonds analyses, and ovals represent final products.

### Pooled data
The main goal of the pooled data analysis is to combine seasonal genetic panels from DEST with new data generated in this study to create a joint pooled-dataset. These data set was filtered to remove samples with high D. simulans contamination as well as for samples with low temporal replication (see paper). The code to do this merger can be found in [DEST](https://dest.bio/), particularly in the mapping pipeline [git](https://github.com/DEST-bio/DEST_freeze1/tree/main/mappingPipeline).

The joint dataset is then used in two major analysis pipelines called "Overwintering" and "GLM". Overwintering consists in a series of multivariate analyses as well as Fst test to assess bottlenecks resulting from boom-and-bust demography and their signal in our data. These data are further filtered according to the effective coverage of the pools.  These results are then used as the basis of a simulation using [SLIM](https://messerlab.org/slim/). The GLM pipeline, on the other hand, uses all the seasonal pooled data to fit model of allele frequency change as a function of weather data obtained from the [NASA POWER](https://power.larc.nasa.gov/) dataset. In this analyses, no effective coverage filter is used because effective coverage is used as a weighting parameter in the models. These models are fit for the Virginia data as well as for other population clusters in DEST. The output of these models are the candidate SNPs for seasonality.

```mermaid
graph LR
A[DEST v1.0] -- Kapun et al --> B[Pools Joints: Code 2.0]
C[VA data, see 0.0 Collections] -- Code from DEST --> B
B --Filter Neff + Overwinter-->Z[Overwinter Set]
Z --Code 3.0/4.0--> 2{PCA time + space}
Z --Code 3.0/4.0--> 1{Fst overwinter}
rr{simulations}
2 --Code 5.0--> rr
1 --Code 5.0--> rr
B --Code 6.0--> E{GLM}
F[NASA POWER: Weather] --Code 6.0--> E
E --Code 6.0--> 22(Best Model VA)
E --Code 6.0--> 33(Best Models EU-W, EU-E, NoA-E)
22 --> 11(Candidate SNPs VA)
33 --> 44(Candidate SNPs Others)
Z --Code 11.0--> gg{Inversion Fst} 
22 --Code 11.0--> gg
```

### Individual data
```mermaid
graph LR
00[Ind. Joint]
AA[VA Individuals] --Code 1.0--> 00
BB[DGRP] --> 00
DD[DPGP] --> 00
CC[Other inds.] --> 00
00 --Code 10.0--> 22{Haplotype Plots}
00 --> 11{Pi, D, Fst}
00 --Code 9.0--> 66{TMRCA}
AA --Code 1.0--> 33{VA LD: Code 8.0}
zz(Candidate SNPs VA  **Pooled dat**) --from Pooled data--> 44
33 --Code 8.0--> 44(Anchor SNPs LDs)
44 --Code 12.0--> qq{Allele Trajectory}
99[Dat from Sanches-Refusta] --Code 12.1--> qq
F[NASA POWER] --weather--> 44
BB --Code 7.0--> pp[In2Lt markers]
pp --> 33
kk{SVM predict Inv}
kk --Code 7.1--> 00
pp --Code 7.1--> kk 
```
### Phenotype analysis: all in Code 13
```mermaid
graph LR
1[226 phenotype sets]
2[DGRP]
3[Geno + Pheno]
4{Inversion GLM effect}
5(Best Model VA)
6{Meta GWAS}
7{Seasonal Trait PCA}
8[Geno+Pheno - Seas. subset]
9{Deficiency Validation}

3 --> 6
5 --> 6
1 --> 3
2 --> 3
3 --> 4
3 --> 8
8 --> 7
5 --> 8
8 --> 9
```

### MSP300 Case Study:  all in Code 14
```mermaid

```

# Files

There are multiple files needed to reproduce our analysis. These files are all publicly available and can be downloaded from the following sites.
