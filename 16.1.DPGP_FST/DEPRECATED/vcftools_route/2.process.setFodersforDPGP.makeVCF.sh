## untar files from DPGP3
tar -xvf dpgp3_sites_vcfs.tar

### make file vector
files=$(ls | grep "vcf.bz2")

for i in ${files[@]}
do
echo $i
mkdir ${i}_fold
mv $i ${i}_fold
cp ./RefGenome.idx.bz2 ./${i}_fold
done