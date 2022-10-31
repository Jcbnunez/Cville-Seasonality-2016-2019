#This is a basic demographic model in msPrime

import numpy as np
import msprime
import IPython
from IPython.display import SVG
import os
import math
import sys
import allel

#setting current working directory
os.chdir("./")

#write vcf with given pre-determined scenario

#set up variables
popA_size = 100_000 #population size individuals are drawn from, held constant
sample_num=1500 #individuals sampled
mu=2.8e-9 #mutation rate
rr=1.25e-7 #recombination rate
seq_len=1e6

#var = argv[1]
var = str('Neutral_100k')

#write demography
demography = msprime.Demography()
demography.add_population(name="A", initial_size=popA_size)

#simulate demographic history
ts = msprime.sim_ancestry(samples=sample_num, demography=demography, sequence_length=seq_len, recombination_rate=rr, model="hudson", ploidy=2)
#simulate random mutations during history
mts = msprime.sim_mutations(ts, rate=mu)

#writes vcf
with open("%s.vcf" % var, "w") as vcf_file:
                      mts.write_vcf(vcf_file)
                      
