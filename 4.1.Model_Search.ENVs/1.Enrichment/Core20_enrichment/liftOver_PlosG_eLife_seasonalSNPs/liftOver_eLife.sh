cat /project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.glm | \
sed '1d' | awk '{print "chr"$1"\t"$2"\t"$2+1"\t"$3"\t"$4"\tdm3_"$1"_"$2}' | tr ' ' '\t' | grep -v "e+" > \
/project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.dm3.bed

~/liftOver \
/project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.dm3.bed \
/scratch/aob2x/dest/dgn/liftoverChains/dm3ToDm6.over.chain \
/project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.dm6.bed \
/project/berglandlab/alan/drosRTEC/mel_all_paired20_2sample_caF_popyear.f_s.unmapped.bed
