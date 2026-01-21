#!/bin/bash
#SBATCH --job-name=Fastq-dump	      # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail  

INDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues"
INPUT_DIR="$INDIR"

ml SRA-Toolkit

cd $INDIR

# Check if the SRR directories exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: File $INPUT_DIR not found. Please ensure it contains SRR directories"
    exit 1
fi

# Loop through each directory in the input directory
for SRR_DIR in "$INPUT_DIR"/*; do
    if [ -d "$SRR_DIR" ]; then
        echo "Processing directory: $SRR_DIR"
        
        # Look for .sra files in the directory
        for SRA_FILE in "$SRR_DIR"/*.sra; do
            if [ -f "$SRA_FILE" ]; then
		    # Extract accession number from the .sra file name (e.g., SRR123456 from SRR123456.sra)
                ACCESSION=$(basename "$SRA_FILE" .sra)

                # Check if FASTQ files already exist for this accession
                if [ -f "$OUTDIR/${ACCESSION}_1.fastq.gz" ] || [ -f "$OUTDIR/${ACCESSION}_2.fastq.gz" ] || [ -f "$OUTDIR/${ACCESSION}.fastq.gz" ]; then
                    echo "FASTQ files already exist for $ACCESSION. Skipping."
                    continue
                fi

                echo "Found SRA file: $SRA_FILE. Processing $ACCESSION..."
                
                # Run fastq-dump on the .sra file
                fastq-dump --outdir "$OUTDIR" --gzip --split-files "$SRA_FILE"
                
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to process $SRA_FILE. Skipping."
                    continue
                fi

                echo "$SRA_FILE processed successfully."
            else
                echo "No .sra files found in $SRR_DIR. Skipping."
            fi
        done
    else
        echo "Skipping $SRR_DIR as it is not a directory."
    fi
done

echo "All directories processed. FASTQ files are stored in $OUTDIR."
