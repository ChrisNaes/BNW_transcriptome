#!/bin/bash
#SBATCH --job-name=Align_Hisat2       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail

#!/bin/bash

# Define directories and reference genome
INDIR="/scratch/cn68176/RNAseq_BNW/Trimmed_reads"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Alignment"
REFDIR="/scratch/cn68176/RNAseq_BNW/Genome" # Adjust depending on the genome needed
REF="vu-2k.fasta" # Adjust depending on the genome needed

# Load necessary modules
ml Miniforge3/24.7.1-0
ml Miniconda3
source activate /home/cn68176/hisat2_conda

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

