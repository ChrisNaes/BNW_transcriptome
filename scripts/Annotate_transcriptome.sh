# Define directories and reference genome
INDIR="Genome"
OUTDIR="Genome/Tissues/Transcriptome"
# SETUP
GENOME_FA="$INDIR/vu-2k.fasta"        # Path to your genome FASTA
THREADS=8                      # Adjust for your machine
GTF_LIST="$INDIR/gtf_list.txt"        # Will be created automatically

ml StringTie

mkdir -p "$INDIR/gtf_out" "$INDIR/merged"
> "$GTF_LIST"

# STEP 1: Assemble transcripts from each sample
#echo "Step 1: Running StringTie per sample..."
for BAM in $OUTDIR/*_sorted.bam; do
    BASENAME=$(basename "$BAM" .bam)
    OUT_GTF="$INDIR/gtf_out/${BASENAME}.gtf"
    echo "  Processing $BAM -> $OUT_GTF"
    stringtie "$BAM" -p $THREADS -o "$OUT_GTF" -l "$BASENAME"
    echo "$OUT_GTF" >> "$GTF_LIST"
done

# STEP 2: Merge all transcript assemblies into unified annotation
echo "Step 2: Merging transcripts..."
stringtie --merge -p $THREADS -o $INDIR/merged/merged.gtf "$GTF_LIST"

# STEP 3: Optional comparison to known annotation (if you have one)
 echo "Step 3: Comparing to known annotation (optional)..."
 gffcompare -r known_reference.gtf -o merged/compare merged/merged.gtf

# STEP 4: Use merged.gtf for downstream quantification
echo "Annotation complete. GTF file is at merged/merged.gtf"


