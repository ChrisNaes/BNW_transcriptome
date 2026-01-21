#!/bin/bash
#SBATCH --job-name=stringtie       # Job name
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
INDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues/Transcriptome1"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Stringtie_wRef"
GTF="/scratch/cn68176/RNAseq_BNW/Genome"

ml StringTie

for FILE in $INDIR/*_sorted.bam; do
       #Extract the base filename
       base_name=$(basename "$FILE" _sorted.bam)
       # assembl the transcriptomes one by one.
       stringtie $INDIR/${base_name}_sorted.bam -o $OUTDIR/${base_name}.gtf -G $GTF/combined.dedup.fixed2.gtf
done
