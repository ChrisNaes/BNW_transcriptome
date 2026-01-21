#!/bin/bash
#SBATCH --job-name=Fastq-dump	      # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail  

INDIR=/scratch/cn68176/RNAseq_H5N1/Raw_reads/Macaque/
OUTDIR=/scratch/cn68176/RNAseq_H5N1/Raw_reads/Macaque/

ml SRA-Toolkit

for file in "$INDIR"/*; do
    # Extract the filename without the extension
    filename=$(basename "$file")

    # Run fastq-dump on all SRR to extract the fastq
    fastq-dump "$filename" -o $OUTDIR/ 
done
