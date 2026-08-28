# Define directories and reference genome
INDIR="Genome/Tissues/Trimmed"
OUTDIR="Genome/Tissues/Transcriptome1"
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

