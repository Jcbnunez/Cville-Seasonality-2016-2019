Some pointers: only .sh files need be run, they call everything else required

Pipeline:

1. poolgens.sh: this just generates the input files for a BayPass run with the same SNPs used in the GLM.
2. matgens.sh: this generates the thinned matrixes used later on. 
These first 2 include some deprecated LOCO features that can be removed
3. wholestand.sh: this generates our actual empirical data, including XtX, BF, and unused contrast statistics
4. simarray.sh: this generates our PODs and resulting BayPass analyses, with all covariables
5. tempmaxarray.sh: this runs BayPass with only tempmax covariable to check for BF convergence.
6. validationsims.sh: this generates PODs and analyses for only tempmax covariable, to allow us to construct a null distribution.
