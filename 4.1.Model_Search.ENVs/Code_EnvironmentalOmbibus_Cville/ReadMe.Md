1. The raw output can be found here: `/project/berglandlab/alan/environmental_ombibus/rawData`. There are 5000 files. Each file has a slice of the genome with ~132 variables_windows x 100 permutations

2. Each variable by window combination has a folder. For example, `/project/berglandlab/alan/environmental_ombibus/temp.ave_4` would be the average temperature using "model 4" which is model F which is 7-15. Right now, inside that folder is a single file, e.g. `temp.ave_4.glmRNP.Rdata` which contains the whole genome for perms 0-100. Rank normalized pvalues are also calculated within each permutation (or within the real data).

3. There is a summary file `/project/berglandlab/alan/environmental_ombibus/bestAIC.Rdata` which contains a datatable. One entry per SNP per permutation (or real). This contains the name of the model with the lowest AIC as well as the AIC values for that model and for the year model.

4. The repository for these scripts is sort of tidy, and see the readme here `Overwintering_18_19/EnvironmentalOmbibus_Cville`



### scripts

1. `Overwintering_18_19/EnvironmentalOmbibus_Cville/getEnvironmentMatrix.cville.R` Generates environmental information matrix
2. `Overwintering_18_19/EnvironmentalOmbibus_Cville/environmentalOmbinus.cville.R` Runs the GLM
