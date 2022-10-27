### artificially create diploids

ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem

module load tabix
module load vcftools

#cp /scratch/yey2sn/Supergene_paper/2.DPGP/SEQs/dpgp3_sequences/DPGP3.2L.merged.vcf ./

vcf_base=/scratch/yey2sn/Supergene_paper/2.DPGP/SEQs/dpgp3_sequences/DPGP3.2L.merged.vcf.gz

##bgzip -d $vcf_base
#head -n 100 DPGP3.2L.merged.vcf > test.vcf
#head -n 100000 DPGP3.2L.merged.vcf > test.vcf

### diploidize

sed 's|\t1|\t1/1|g' DPGP3.2L.merged.vcf | \
sed 's|\t0|\t0/0|g' | \
sed 's|\t2|\t2/2|g' \
> DPGP3.2L.merged.fixGT.vcf


####
KAR=DPGP3
#File containing populations
std_samps=std.dpgp3.txt
inv_samps=inv.dpgp3.txt

window2=500000
step2=100000

window3=10000
step3=5000

window4=5000
step4=1000


#######

vcftools \
--vcf DPGP3.2L.merged.fixGT.vcf \
--weir-fst-pop $std_samps \
--weir-fst-pop $inv_samps \
--fst-window-size $window3 \
--fst-window-step $step3 \
--chr chr2L \
--out FST.W_$window3.S_$step3.$KAR
