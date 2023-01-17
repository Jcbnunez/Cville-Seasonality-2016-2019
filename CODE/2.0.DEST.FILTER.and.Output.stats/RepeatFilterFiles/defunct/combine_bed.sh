module load gcc/9.2.0 bedtools/2.29.2

cd ~/

cat \
Overwintering_18_19/RepeatFilterFiles/InterruptedRepeats.bed \
Overwintering_18_19/RepeatFilterFiles/MicroSats.bed \
Overwintering_18_19/RepeatFilterFiles/RepeatMasker.bed \
Overwintering_18_19/RepeatFilterFiles/SimpleRepeats.bed \
Overwintering_18_19/RepeatFilterFiles/WM_SDust.bed |
grep -v "track name" | cut -f1-4 > Overwintering_18_19/RepeatFilterFiles/combined_repeats.bed

bedtools sort \
-i Overwintering_18_19/RepeatFilterFiles/combined_repeats.bed | \
bedtools merge > \
Overwintering_18_19/RepeatFilterFiles/combined_repeats.sort.merge.bed
