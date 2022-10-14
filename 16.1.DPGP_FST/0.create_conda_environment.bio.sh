module load anaconda

conda create \
-n bio_space \
python=3.8 

source activate bio_space

pip install bio --upgrade
