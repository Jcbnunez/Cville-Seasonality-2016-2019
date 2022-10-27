# A cosmopolitan inversion drives seasonal adaptation in overwintering Drosophila

This is the git repo for the paper "A cosmopolitan inversion drives seasonal adaptation in overwintering Drosophila". This repo contains all the code and data needed to replicate our analyses.

# Citation
Please cite our paper as TBD

# Files

There are multiple files needed to reproduce our analysis. These files are all publicly available and can be downloaded from the following sites.

## Analysis scheme flowchart

### Pooled data
```mermaid
graph LR
A[DEST v1.0] -- pool-seq --> B[Pools Joints]
C[VA data] -- pool-seq --> B
B --Filter Neff + Overwinter-->Z[Overwinter Set]
Z --> 2{PCA time + space}
Z --> 1{Fst overwinter}
B --> E{GLM}
F[NASA POWER] --weather--> E
E --> 22(Best Model VA)
E --> 33(Best Models EU-W, EU-E, NoA-E)
22 --> 11(Candidate SNPs VA)
33 --> 44(Candidate SNPs Others)
```

### Individual data
```mermaid
graph LR
00[Individuals Joint]
AA[VA Individuals] --> 00
BB[DGRP] --> 00
CC[Other ind dat. see below] --> 00
00 --> 22{Haplotype Plots}
00 --> 11{Pi, D, Fst}
AA --> 33{VA LD}
zz(Candidate SNPs VA - from pools) --> 33
33 --> 44(Anchor SNPs LDs)
44 --> qq{Allele Trajectory}
F[NASA POWER] --weather--> 44
```

### Joint analyses 
```mermaid
graph LR
00[Individuals Joint]

```
