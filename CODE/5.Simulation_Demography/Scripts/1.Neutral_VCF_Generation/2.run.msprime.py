# This is a basic demographic model in msPrime
# Generate neutral VCF for SLIM

# Load modules
import numpy as np
import msprime
import IPython
from IPython.display import SVG
import os
import math
import sys
import allel

# Working directory
os.chdir("/project/berglandlab/connor/slim_bottleneck/")

# Model parameters
popA_size = int(sys.argv[1]) # Constant population size
sample_num = 1500 # End point individuals sampled
mu = 2.8e-9 # Mutation rate
rr = 1.25e-7 # Recombination rate
seq_len = 1e6 # Length of chromosome

# Array name
var = str(sys.argv[2])

# Prewritten demographic model
demography = msprime.Demography()
demography.add_population(name="A", initial_size=popA_size)

# Simulate demographic history
ts = msprime.sim_ancestry(samples=sample_num,
                          demography=demography,
                          sequence_length=seq_len,
                          recombination_rate=rr,
                          model="hudson",
                          ploidy=2)

# Simulate mutations
mts = msprime.sim_mutations(ts, rate=mu)

# Outputs VCF
with open("%s.vcf" % var, "w") as vcf_file:
                      mts.write_vcf(vcf_file)
