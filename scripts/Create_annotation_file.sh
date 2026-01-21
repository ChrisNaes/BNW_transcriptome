#!/bin/bash
#SBATCH --job-name=liftoff			 # Job name
#SBATCH --partition=batch        		 # Partition (queue) name
#SBATCH --ntasks=100              		 # Run on a single CPU
#SBATCH --mem=20gb                		 # Job memory request
#SBATCH --time=100:00:100         		 # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out     # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err      # Standard error log
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu              # Where to send mail

INDIR="/scratch/cn68176/RNAseq_BNW/Genome"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Stringtie"
REF="vu-2k.fna" # Adjust depending on the genome needed

ml Miniforge3

source activate liftoff 

liftoff -g $INDIR/GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.gtf \
	-o $INDIR/liftoff_vu2k.gtf \
	-p 8 -a 0.9 -s 0.9 \
	$INDIR/vu-2k.fasta $INDIR/GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.fna
