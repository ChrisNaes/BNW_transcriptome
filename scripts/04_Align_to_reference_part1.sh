#!/bin/bash
#SBATCH --job-name=Align_Hisat2       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail

INDIR="/scratch/cn68176/RNAseq_H5N1/Trimmed_reads/RiChicken"
OUTDIR="/scratch/cn68176/RNAseq_H5N1/Alignment/RiChicken"
REFDIR="/scratch/cn68176/RNAseq_H5N1/Genomes/Chicken_Gallus_gallus" # Adjust depending on the genome needed
REF="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna" # Adjust depending on the genome needed

ml Miniforge3/24.7.1-0
ml Miniconda3

source activate /home/cn68176/hisat2_conda

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

