
#load module
module load plink


#load input data
vcf=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/ForPlinkMarkersOnly.CM.2L.recode.vcf.gz

ld_snp=2L_12899672_SNP
temp_snp=2L_748285_SNP

runid=`uuidgen`
echo $runid

#run plink
plink --vcf $vcf \
--allow-extra-chr \
--double-id \
--ld $ld_snp $temp_snp \
--out tmp.$ld_snp.$temp_snp

#clean up
rm tmp.$ld_snp.$temp_snp.nosex

#slice out the valuable information
cat tmp.$ld_snp.$temp_snp.log | \
sed '1,33d' | \
head -n 4 | \
sed -E 's/\ +/\ /g' | \
sed -E 's/\ /\t/g' | \
sed -n '1p;3p' > tmp.$ld_snp.$temp_snp.$temp_snp.flt.txt

rm tmp.$ld_snp.$temp_snp.log

paste -d "\t" tmp.$ld_snp.$temp_snp.flt.txt <(printf %"s\n" "present" "absent") \
> tmp.$ld_snp.$temp_snp.flt.hdr.txt

rm tmp.$ld_snp.$temp_snp.flt.txt

sed -i "s/^/${ld_snp}vs${temp_snp}\t/" tmp.$ld_snp.$temp_snp.flt.hdr.txt



