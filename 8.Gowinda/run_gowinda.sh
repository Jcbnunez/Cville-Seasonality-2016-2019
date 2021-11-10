Gowinda=/home/yey2sn/software/Gowinda/Gowinda-1.12.jar
Universe=glm_universe.2L.inversion.txt
Outliers=glm_p01_targets.p001.2L.inversion.txt
Assocs=/project/berglandlab/Dmel_genomic_resources/Annotations/funcassociate_go_associations_CGstyle.txt
Annots=/project/berglandlab/Dmel_genomic_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_CGStyle.gtf


###
# Make GTF
#gff=/project/berglandlab/Dmel_genomic_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff

#Gff2Gtf=/home/yey2sn/software/Gowinda/Gff2Gtf.py
#python $Gff2Gtf --input $gff > GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf

java -Xmx35g \
-jar $Gowinda \
--snp-file $Universe \
--candidate-snp-file $Outliers \
--gene-set-file $Assocs \
--annotation-file $Annots \
--simulations 100000 \
--min-significance 1 \
--gene-definition updownstream2000 \
--threads 10 \
--output-file GLMtemp_results_gene_SNP_2kups.2LtInv.txt \
--mode snp \
--min-genes 1



