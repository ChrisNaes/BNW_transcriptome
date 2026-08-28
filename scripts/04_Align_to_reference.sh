INDIR="Trimmed_reads"
OUTDIR="Alignment"
REFDIR="Genomes" # Adjust depending on the genome needed
REF="vu-2k.fna" # Adjust depending on the genome needed

ml Miniforge3/24.7.1-0
ml Miniconda3

source activate hisat2_conda

# Navigate to the directory where the reference genome files are located
cd $REFDIR

# Ensure the reference genome file is unzipped
if [ -f "${REF}.gz" ]; then
    gunzip -c "${REF}.gz" > "${REF}"
fi

# Ensure the reference genome is indexed
if [ ! -f "${REF}.1.ht2" ]; then
    hisat2-build -f "${REF}" "${REF}"
fi

------------------------------------------------

# Set number of threads (adjust as needed)
THREADS=8  # Change this based on available CPU resources

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Function to run HiSat2 alignment in parallel
align_reads() {
    sample=$(basename "$1" _R1_001_trimmed.fq.gz)  # Extract sample name
    hisat2 -p "$THREADS" -x "${REFDIR}/${REF}" -U "${INDIR}/${sample}_R1_001_trimmed.fq.gz" -S "${OUTDIR}/${sample}.sam"
}

export -f align_reads
export INDIR OUTDIR REFDIR REF THREADS

# Run HiSat2 alignment in parallel using GNU parallel

ls ${INDIR}/*_R1_001_trimmed.fq.gz | sed 's/_R1_001_trimmed.fq.gz//' | xargs -I{} -P 4 bash -c 'align_reads "$@"' _ {}

# Explanation:
# - `ls ${INDIR}/*_1.fq` lists all forward read files
# - `sed 's/_1.fq//'` removes `_1.fq` to get the base sample name
# - `xargs -I{} -P 4 bash -c 'align_reads "$@"' _ {}` runs `align_reads` in parallel with 4 jobs at a time (`-P 4`)
# - Adjust `-P 4` to optimize CPU usage
