#!/bin/bash
#SBATCH --ntasks-per-node=15
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez
#SBATCH --array=1-5


## link main files related to the genome of D.mel
Gowinda=/home/yey2sn/software/Gowinda/Gowinda-1.12.jar
Assocs=/project/berglandlab/Dmel_genomic_resources/Annotations/funcassociate_go_associations_CGstyle.txt
Annots=/project/berglandlab/Dmel_genomic_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_CGStyle.gtf

# Different universes of SNPs
#Universe=glm_universe.2L.inversion.txt
#Outliers=glm_p01_targets.p001.2L.inversion.txt

### Optional Step to make the GTF file... this often is needed only once
# Make GTF
#gff=/project/berglandlab/Dmel_genomic_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff

#Gff2Gtf=/home/yey2sn/software/Gowinda/Gff2Gtf.py
#python $Gff2Gtf --input $gff > GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf


guide=guide_file.txt

#extract pop name from guide file
i=$( sed -n ${SLURM_ARRAY_TASK_ID}p $guide |  awk '{print $1}' )
k=$( sed -n ${SLURM_ARRAY_TASK_ID}p $guide |  awk '{print $2}' )

echo $i
echo $k
	
	java -Xmx35g \
	-jar $Gowinda \
	--snp-file ${i}_noTresh_glm_universe.txt \
	--candidate-snp-file ${i}_${k}_glm_outliers.txt \
	--gene-set-file $Assocs \
	--annotation-file $Annots \
	--simulations 100000 \
	--min-significance 1 \
	--gene-definition updownstream2000 \
	--threads 10 \
	--output-file GOWINDA_out_${i}_${k}_SNPmode_100kperms.txt \
	--mode snp \
	--min-genes 1 \
	--detailed-log

echo "done"
