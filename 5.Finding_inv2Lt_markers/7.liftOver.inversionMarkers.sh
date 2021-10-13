# Lift over code to transform dm3 inversion markers to dm6


#load the liftOver program
liftOver=/project/berglandlab/alan/liftOver

#define variables
in_file=inv_2l_dm3_markers_asPos.txt
type=pos

chain_file=/project/berglandlab/Dmel_genomic_resources/liftOver_files/dm3ToDm6.over.chain

#define the outputfiles
out_file=inversion_2lt_makers.Dm6.$type.bed
out_unmapped=inversion_2lt_makers.Unmapped.Dm6.$type.bed


################
#Run lift over
$liftOver -positions \
$in_file \
$chain_file \
$out_file \
$out_unmapped

head $out_file
wc -l $out_unmapped
