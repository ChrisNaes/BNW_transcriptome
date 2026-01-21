#!/bin/bash
#SBATCH --job-name=annotate       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail

# Define directories and reference genome
INDIR="/scratch/cn68176/RNAseq_BNW/Genome"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues/Transcriptome"
# SETUP
GENOME_FA="$INDIR/vu-2k.fasta"        # Path to your genome FASTA
THREADS=8                      # Adjust for your machine
GTF_LIST="$INDIR/gtf_list.txt"        # Will be created automatically

ml StringTie

mkdir -p "$INDIR/gtf_out" "$INDIR/merged"
> "$GTF_LIST"

# STEP 1: Assemble transcripts from each sample
#echo "Step 1: Running StringTie per sample..."
#for BAM in $OUTDIR/*_sorted.bam; do
#    BASENAME=$(basename "$BAM" .bam)
#    OUT_GTF="$INDIR/gtf_out/${BASENAME}.gtf"
#    echo "  Processing $BAM -> $OUT_GTF"
#    stringtie "$BAM" -p $THREADS -o "$OUT_GTF" -l "$BASENAME"
#    echo "$OUT_GTF" >> "$GTF_LIST"
#done

# STEP 2: Merge all transcript assemblies into unified annotation
echo "Step 2: Merging transcripts..."
stringtie --merge -p $THREADS -o $INDIR/merged/merged.gtf "$GTF_LIST"

# STEP 3: Optional comparison to known annotation (if you have one)
# echo "Step 3: Comparing to known annotation (optional)..."
# gffcompare -r known_reference.gtf -o merged/compare merged/merged.gtf

# STEP 4: Use merged.gtf for downstream quantification
#echo "Annotation complete. GTF file is at merged/merged.gtf"


