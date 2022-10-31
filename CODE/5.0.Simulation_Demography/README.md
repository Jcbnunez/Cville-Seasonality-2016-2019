# 5.0.Simulation_Demography

## 0.Simulate_Neutral_Demography
Coalescent Simulation work

### 0.install.msprime.sh
This code sets up msprime in the environment (conda)
### 1.run.msprime.py
Generate a netral demography to serve as burn-in for the forward simulation.


## 1.Simulate_bottlenecks
Forward Simulation work

### 0.model_parameter_generator.R
Aux code to generate parameters for the simulations

### 1.bottleneck_discrete.slim
Core SLIM code that simulates the boom-bust demography

### 2.bottleneck_discrete_reps_forloop_readvcf.sh
Launcher of 1.bottleneck_discrete.slim

### 3.slim_vcf_output_readvcf_pairwisefst_pca.R
Parse the output of the simulations

### 4.analysis_slim_fst_new3.R
First round of analysis for the output data

### 5.analysis.abc.R
Preliminary runs of ABC on the data (there a folder below with the formal analysis)

## 2.ABC_est
ABC test

### 1.Parse.data.for.ABC.R
Prepare the data for ABC

### 2.Launch.Parser.sh
Launch 1.Parse.data.for.ABC.R

### 3. join.datasets.DoABC.R
Joint real and simulated data for ABC

### 4.calculate.harmonic.ne.R
Calculate harmonic Ne (estimated)  in the population
