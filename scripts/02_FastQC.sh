#!/bin/bash
#SBATCH --job-name=Quality_check      # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail  

INDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues"
OUTDIR="/scratch/cn68176/RNAseq_BNW/QC"

cd $INDIR

ml FastQC/0.11.9-Java-11
ml MultiQC/1.14-foss-2022a

# Loop over all fq.gz files in the input directory
for file in "$INDIR"/*.fastq.gz; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fastq.gz)

    # Run FastQC on all trimmed reads to make sure everything looks okay
    fastqc "$filename".fastq.gz -o $OUTDIR/ 
done

# Go into the directory with the QC log files
cd $OUTDIR

# Run MultiQC to summarise the output from FastQC. 
# -d command tells multiqc to look in all potential folders to find log files. 
multiqc .  
