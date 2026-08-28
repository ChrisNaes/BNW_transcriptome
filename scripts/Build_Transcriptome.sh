# Define directories and reference genome
INDIR="Genome/Tissues/Trimmed"
OUTDIR="Genome/Tissues/Transcriptome"
REFDIR="Genome" # Adjust depending on the genome needed
REF="vu-2k.fasta" # Adjust depending on the genome needed

# Load necessary modules
ml Miniforge3/24.7.1-0
ml Miniconda3
source activate hisat2_conda

hisat2-build $REFDIR/$REF $REFDIR/$REF

# Set number of threads (adjust as needed)
THREADS=8  # Change this based on available CPU resources

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Function to run HiSat2 alignment in parallel
align_reads() {
    sample=$(basename "$1" _R1_val_1.fq.gz)  # Extract sample name
    hisat2 --phred33 --rna-strandness RF -p "$THREADS" --novel-splicesite-outfile "$OUTDIR/${sample}_splicesite.txt" -S "$OUTDIR/${sample}_accepted_hits.sam" -x "${REFDIR}/${REF}" -1 "${INDIR}/${sample}_R1_val_1.fq.gz" -2 "${INDIR}/${sample}_R2_val_2.fq.gz" -S "${OUTDIR}/${sample}.sam"
}

export -f align_reads
export INDIR OUTDIR REFDIR REF THREADS

# Run HiSat2 alignment in parallel using GNU parallel
ls ${INDIR}/*_R1_val_1.fq.gz | sed 's/_R1_val_1.fq.gz//' | xargs -I{} -P 4 bash -c 'align_reads "$@"' _ {}

-----------------------------------------------------------------

ml SAMtools

## We first use samtools to create sorted and indexed BAM files for our samples
for FILE in $OUTDIR/*.sam; do
       #Extract the base filename
       base_name=$(basename "$FILE" .sam)

       # Convert SAM to BAM files
       samtools view -bS "$OUTDIR/${base_name}.sam" > "$OUTDIR/${base_name}.bam"

done

cd $OUTDIR
# Now merge the tissue specific bam files into one per tissue
samtools merge -o Liver.bam Liver_S19_L003.bam Liver_S25_L004.bam
samtools merge -o Lung.bam Lung_S22_L003.bam Lung_S28_L004.bam
samtools merge -o Lymph-Node.bam Lymph-Node_S18_L003.bam Lymph-Node_S24_L004.bam
samtools merge -o Skeletal-Muscle.bam Skeletal-Muscle_S20_L003.bam Skeletal-Muscle_S26_L004.bam 
samtools merge -o Skin.bam Skin_S17_L003.bam Skin_S23_L004.bam
samtools merge -o Thyroid.bam Thyroid_S21_L003.bam Thyroid_S27_L004.bam 

for FILE in $OUTDIR/*.bam; do
#       #Extract the base filename
       base_name=$(basename "$FILE" .bam)
#       # Sort BAM files
       samtools sort -@ 8 "$OUTDIR/${base_name}.bam" -o "$OUTDIR/${base_name}_sorted.bam"

#       # Index the sorted BAM files
       samtools index "$OUTDIR/${base_name}_sorted.bam"

#       # Print the information for the sample that completes then remove the temporary files
       echo ${base_name}.bam
done
