cat /project/berglandlab/alan/drosRTEC/seas2014.delim | \
awk '{print "chr"$1"\t"$2"\t"$2+1"\t"$3}' | tr ' ' '\t' | grep -v "e+" > \
/project/berglandlab/alan/drosRTEC/seas2014.dm3.bed

~/liftOver \
/project/berglandlab/alan/drosRTEC/seas2014.dm3.bed \
/scratch/aob2x/dest/dgn/liftoverChains/dm3ToDm6.over.chain \
/project/berglandlab/alan/drosRTEC/seas2014.dm6.bed \
/project/berglandlab/alan/drosRTEC/seas2014.unmapped.bed
