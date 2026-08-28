
# Define directories and reference genome
INDIR="Genome/Tissues/Transcriptome"
OUTDIR="Genome/Stringtie_wRef"
GTF="Genome"

ml StringTie

for FILE in $INDIR/*_sorted.bam; do
       #Extract the base filename
       base_name=$(basename "$FILE" _sorted.bam)
       # assembl the transcriptomes one by one.
       stringtie $INDIR/${base_name}_sorted.bam -o $OUTDIR/${base_name}.gtf -G $GTF/combined.dedup.fixed.gtf
done
