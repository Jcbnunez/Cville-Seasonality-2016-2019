# 1.0.Map_Individual_DNAseqs

## Files 
Need to download reads from the SRA

## Make pair files
ls  /files/fastq/ | grep "HL5KFALXX" > HL5KFALXX.file_names.txt

## Make a file with all the file pairs from the experiment

head ../all_mapping_pairs.txt

## Check the data
Notice that not all the samples of the SRA were used due to quality reasons. Use only the files described in Table S1.

see ../Check_mapping_files.R

### USAGE NOTES -- Using a SLURM manager

```
echo "have you checked if the folder for SLURMout already built with mkdir?"
if [[ -d "SLURMout" ]]
then
	echo "Working merged_reads folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir ./SLURMout
	date
fi
```

### Step 1: Merge and QC
Merge paired end reads into large reads and do QC
```
sbatch --array=1-$( wc -l < ./Alyssa_ind_reads_guideFile.txt ) \
 --output ./SLURMout/Merge_and_QC_reads.%A_%a.out \
 --error ./SLURMout/Merge_and_QC_reads.%A_%a.err \
./1.Merge_and_QC_reads.2016samps.sh
```
### Step 1.1: Run Merge fast QC

```
module load singularity
singularity pull --name multiqc.sif shub://cory-weller/singularity:multiqc
singularity run /home/yey2sn/software/multiqc.sif ./fastqc_merged
```
### Step 2: Trim and Map
Trim read edges and Map to Drosophila hologenome using `bwa`

```
sbatch --array=1-$( wc -l < ./Alyssa_ind_reads_guideFile.txt ) \
 --output ./SLURMout/Trim_and_Map.%A_%a.out \
 --error ./SLURMout/Trim_and_Map.%A_%a.err \
./2.Trim_and_Map.2016samps.sh
```

### Step 3 Merge bams across lanes
Merge bam files from different mapping strategies
```
sbatch --array=1-$( wc -l < ./Alyssa_ind_logic_comparison_inds.txt ) \
 --output ./SLURMout/MergeLanes.%A_%a.out \
 --error ./SLURMout/MergeLanes.%A_%a.err \
./3.MergeBams.2016samps.sh
```

### Step 3.1 -- Run Qualimap
Run qualimap on the bam files
```
module load qualimap
qualimap  multi-bamqc  \
-outdir ./multi_bamQC \
-d ./multi_qualimap.txt
```
### Step 4 -- Run Haplocaller
Haplotype samples with GATK
```
sbatch --array=1-$( wc -l < ./Alyssa_ind_logic_comparison_inds.txt ) \
 --output ./SLURMout/HaploCaller.%A_%a.out \
 --error ./SLURMout/HaploCaller.%A_%a.err \
./4.HaploCaller.2016samps.sh
```

### Step 5 -- Build GenoDB
Build genotyoe database
```
sbatch --array=1-$( wc -l < ./Intervals_Dmel.txt ) \
 --output ./SLURMout/GenoDB.%A_%a.out \
 --error ./SLURMout/GenoDB.%A_%a.err \
./5.Build_GenomeDB.sh
```

### Step 6 -- Genotype Files
Genotype samples
```
sbatch --array=1-$( wc -l < ./Intervals_Dmel.txt ) \
 --output ./SLURMout/Genotype.%A_%a.out \
 --error ./SLURMout/Genotype.%A_%a.err \
./6.GenotypeGVCFs.sh
```

### Step 7 -- Recalibrate VCF
Recalibrate VCF with DGRP
```
sbatch --array=1-$( wc -l < ./Intervals_Dmel.txt ) \
 --output ./SLURMout/Recalibrate.%A_%a.out \
 --error ./SLURMout/Recalibrate.%A_%a.err \
./7.RecalibrateVCF.sh
```
### Step 8 -- merge VCF
Merge VCFs from different chromosomes
```
sbatch \
 --output ./SLURMout/merge.%A_%a.out \
 --error ./SLURMout/merge.%A_%a.err \
./8.MergeVCF.sh
```
### Step 9 -- index all bam files before phasing
```
sbatch --array=1-$( wc -l < ./samples_to_index.txt ) \
 --output ./SLURMout/index.%A_%a.out \
 --error ./SLURMout/index.%A_%a.err \
./9.index_bams.sh
```
### Step 10 -- Phase with Whatshapp
Phase using read information
```
sbatch --array=1-$( wc -l < ./samples_to_phase_chr.txt ) \
 --output ./SLURMout/Whatshapp.%A_%a.out \
 --error ./SLURMout/Whatshapp.%A_%a.err \
./10.Phase_whatshap.sh
```

### Step 11 -- parse Whatshapp output
```
sbatch --array=1-$( wc -l < ./Intervals_for_phasing.txt ) \
 --output ./SLURMout/parse_Whatshapp.%A_%a.out \
 --error ./SLURMout/parse_Whatshapp.%A_%a.err \
./11.Merge_Phase_whatshap.sh
```
### Step 12 -- remove single repeats
```
sbatch --array=1-$( wc -l < ./Intervals_for_phasing.txt ) \
 --output ./SLURMout/removeSR.%A_%a.out \
 --error ./SLURMout/removeSR.%A_%a.err \
./12.Remove_Simple_repeats.sh
```

### Step 13 -- Filter individual outliers using preliminary PCA
```
sbatch \
 --output ./SLURMout/FilterPCA.%A_%a.out \
 --error ./SLURMout/FilterPCA.%A_%a.err \
./13.Filter_VCF_by_PCA.sh
```

### Step 14 -- Phase w/ ShapeIt 
```
sbatch --array=1-$( wc -l < ./Intervals_for_phasing.txt ) \
 --output ./SLURMout/parse_shapeit.%A_%a.out \
 --error ./SLURMout/parse_shapeit.%A_%a.err \
./14.phase_w.shapeit4.sh
```

### Step 15 -- merge BCFs and make VCF
```
sbatch \
 --output ./SLURMout/mergePhased.%A_%a.out \
 --error ./SLURMout/mergePhased.%A_%a.err \
./15.merge_phased_bcfs.sh
```

## Guide Files
Metadata used in this pipeline