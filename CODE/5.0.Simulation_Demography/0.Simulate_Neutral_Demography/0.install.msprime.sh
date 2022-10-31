# Install msprime
# 

module load anaconda/2020.11-py3.8

conda create \
-n msprime_env \
python=3.8 \
scikit-allel \
ipykernel \
-c conda-forge

source activate msprime_env

conda install msprime -c conda-forge