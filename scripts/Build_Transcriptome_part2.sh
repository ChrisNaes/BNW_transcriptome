#!/bin/bash
#SBATCH --job-name=samtools       # Job name
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
INDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues/Trimmed"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues/Transcriptome"
REFDIR="/scratch/cn68176/RNAseq_BNW/Genome" # Adjust depending on the genome needed
REF="vu-2k.fasta" # Adjust depending on the genome needed

ml SAMtools


## We first use samtools to create sorted and indexed BAM files for our samples
#for FILE in $OUTDIR/*.sam; do
       #Extract the base filename
#       base_name=$(basename "$FILE" .sam)

       # Convert SAM to BAM files
#       samtools view -bS "$OUTDIR/${base_name}.sam" > "$OUTDIR/${base_name}.bam"

#done
#cd $OUTDIR
# Now merge the tissue specific bam files into one per tissue
#samtools merge -o Liver.bam Liver_S19_L003.bam Liver_S25_L004.bam
#samtools merge -o Lung.bam Lung_S22_L003.bam Lung_S28_L004.bam
#samtools merge -o Lymph-Node.bam Lymph-Node_S18_L003.bam Lymph-Node_S24_L004.bam
#samtools merge -o Skeletal-Muscle.bam Skeletal-Muscle_S20_L003.bam Skeletal-Muscle_S26_L004.bam 
#samtools merge -o Skin.bam Skin_S17_L003.bam Skin_S23_L004.bam
#samtools merge -o Thyroid.bam Thyroid_S21_L003.bam Thyroid_S27_L004.bam 

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

