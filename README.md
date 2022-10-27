# A cosmopolitan inversion drives seasonal adaptation in overwintering Drosophila

This is the git repo for the paper "A cosmopolitan inversion drives seasonal adaptation in overwintering Drosophila". This repo contains all the code and data needed to replicate our analyses.

# Citation
Please cite our paper as TBD

## Analysis scheme flowchart

### Pooled data
```mermaid
graph LR
A[DEST v1.0] -- Kapun et al --> B[Pools Joints: Code 2.0]
C[VA data] -- Code from DEST --> B
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
tt[Collections] --0.0.Collections--> C
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
### Phenotype analysis
```mermaid
graph LR
```
