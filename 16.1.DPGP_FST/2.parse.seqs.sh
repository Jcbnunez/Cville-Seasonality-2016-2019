### prepare files from DPGP
ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

### after downloading the global tar file; unpack each chromosome arm
mkdir raw_seqs
cp /scratch/yey2sn/Supergene_paper/2.DPGP/SEQs/dpgp3_sequences/tars/dpgp3_Chr2L.tar ./raw_seqs
tar -xvf ./raw_seqs/dpgp3_Chr2L.tar 
 
### create vector of files  
cd ./raw_seqs/
files=$(ls * | grep ".seq")
cd ..

### create folder of outputs
mkdir fasta_outs

#### run loop
for i in ${files[@]}
do

echo $i
echo ">"${i} | sed "s/_Chr2L.seq//"  > header.$i.txt
sed "$ s/$/\n/" ./raw_seqs/$i > $i.EOL.fa
cat header.$i.txt $i.EOL.fa  > fasta_outs/$i.named.fasta

rm header.$i.txt
rm $i.EOL.fa

head -c 30 fasta_outs/$i.named.fasta

done
