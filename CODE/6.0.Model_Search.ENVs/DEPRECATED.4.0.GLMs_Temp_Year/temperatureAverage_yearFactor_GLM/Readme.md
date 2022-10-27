# sripts to run temperature GLMs

## 1.X: various implementations of the thermal GLM.
### 1.1: compare `af~tempAve` vs `af~1`, and `af~year_factor` vs `af~1`
### 1.2: compare `af~tempAve+year_factor` vs `af~year_factor`, and `af~year_factor` vs `af~1`
### 1.3: CURRENT version: compare `af~tempAve+year_factor` vs `af~year_factor`, and `af~year_factor` vs `af~1` using quasibinomial error structure. Incudes scripts to scatter jobs, collect jobs, do rank_normalization. Output data is found to `/project/berglandlab/thermal_glm_dest/processedGLM`. The number refers to the permutation; 0=real.

## 2: Inversion analysis. This uses the data step 1. This is becoming defunct?

## 3: Window analysis. Generates various window based statistics. Also makes a nice inversion type plot. Output data can be found in `/project/berglandlab/thermal_glm_dest/windowAnalysis`. 

## 4: kind of defunct
`
