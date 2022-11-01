# 5.0.drosophila_bottleneck_simulations
* 10/31/2022
* Started by Connor S. Murray

Drosophila bottleneck simulations

Genetic differentiation of a neutrally evolving population after consecutive winter bottlenecks.

## Scripts
# 1) Generate the burn-in VCF

This will generate the VCF for a metapopulation to have a theta pi ~ 0.01 (wild Drosophila metapopulation genetic diversity).

Uses msprime/python

# 2) Run the SLiM scripts and output statistics

Independently replicates the simulated universe over 50 generations with 2 bottleneck events.

Uses SLiM/R

# 3) Merge parsed data

Collect data from each simulation run and concatenates it for analysis step.

Uses bash/R

# 4) Analyze, ABC, and plot data

Approximate bayesian computation (ABC) and plotting of results.

Uses R
