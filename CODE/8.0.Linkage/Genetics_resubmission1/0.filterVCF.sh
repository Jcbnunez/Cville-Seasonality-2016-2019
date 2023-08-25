## some filtering
vcftools=/gpfs1/home/j/c/jcnunez/software/vcftools/vcftools
spack load bcftools@1.14

vcf=/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

pop1=/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Single_Individuals/CM.STD.txt
pop2=/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Single_Individuals/CM.INV.txt
pop3=/netfiles/nunezlab/Drosophila_resources/Nunez_et_al_Supergene_paper/Single_Individuals/CM.HET.txt

## calculate D
$vcftools --gzvcf $vcf \
--remove-indels \
--max-alleles 2 \
--max-missing 0.85 \
--maf 0.01 \
--chr 2L \
--thin 10 \
--keep $pop1 \
--keep $pop2 \
--recode --recode-INFO-all \
--out CM.2L.maf05.mis85.thin100.inv.std

bgzip CM.2L.maf05.mis85.thin100.inv.std.recode.vcf
tabix CM.2L.maf05.mis85.thin100.inv.std.recode.vcf.gz
