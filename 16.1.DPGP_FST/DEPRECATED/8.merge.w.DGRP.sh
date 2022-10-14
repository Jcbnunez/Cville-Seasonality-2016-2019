module load bcftools
module load tabix

dgrp="/project/berglandlab/DGRP_freeze2_vcf/DGRP2_freeze2.2L.vcf.gz"
dpgp="./DPGP3.2L.merged.vcf.gz"

bcftools merge \
$dgrp \
$dpgp \
> DPGP3.DGRP.2L.merged.vcf

bgzip DPGP3.DGRP.2L.merged.vcf
tabix DPGP3.DGRP.2L.merged.vcf.gz

