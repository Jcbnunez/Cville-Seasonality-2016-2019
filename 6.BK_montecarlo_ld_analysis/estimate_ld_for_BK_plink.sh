
#load module
module load plink

#load information
focal_type="inv_focus"

#load list of SNPs to iterate 
iterator_snps_list=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/temperature_snps_ids_p1.txt

#load input data
vcf=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_temp_snps/ForPlinkMarkersOnly.CM.2L.recode.vcf.gz
dependent_snps_list=/scratch/yey2sn/Overwintering_ms/old/Fig6_analyze_dependent_snps/temperature_snps_ids_p1.txt

#set a fix marker SNP
focal_snp=2L_12899672_SNP

#loop over temperature outlier snps
mkdir simulations_output

while read dependent_snp
#dependent_snp=2L_748285_SNP
do # <---------- do here

#run plink
plink --vcf $vcf \
--allow-extra-chr \
--double-id \
--ld $focal_snp $dependent_snp \
--out tmp.$focal_snp.$dependent_snp

#clean up
rm tmp.$focal_snp.$dependent_snp.nosex

#slice out the valuable information
cat tmp.$focal_snp.$dependent_snp.log | \
sed '1,33d' | \
head -n 4 | \
sed -E 's/\ +/\ /g' | \
sed -E 's/\ /\t/g' | \
sed -n '1p;2p' > tmp.$focal_snp.$dependent_snp.flt.txt

rm tmp.$focal_snp.$dependent_snp.log

paste -d "\t" tmp.$focal_snp.$dependent_snp.flt.txt <(printf %"s\n" "present" "absent") \
> tmp.$focal_snp.$dependent_snp.flt.hdr.txt

rm tmp.$focal_snp.$dependent_snp.flt.txt

sed -i "s/^/${focal_snp}vs${dependent_snp}/" tmp.$focal_snp.$dependent_snp.flt.hdr.txt
sed -i "s/^/$focal_type\t/" tmp.$focal_snp.$dependent_snp.flt.hdr.txt

cat tmp.$focal_snp.$dependent_snp.flt.hdr.txt >> simulations_output/$focal_snp.marker.BK.txt
rm tmp.$focal_snp.$dependent_snp.flt.hdr.txt

done < $iterator_snps_list # <---------- done


