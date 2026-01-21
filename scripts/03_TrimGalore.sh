#!/bin/bash
#SBATCH --job-name=TrimGalore		  # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=20gb                    # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out            # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err             # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   # Where to send mail	

INDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues"
OUTDIR="/scratch/cn68176/RNAseq_BNW/Genome/Tissues/Trimmed"

ml Trim_Galore/0.6.7-GCCcore-11.2.0 
ml cutadapt/4.5-GCCcore-11.3.0

# --output_dir specifies the output directory that the trimmed files will be located in. 
# --quality 20 removes low quality ends from reads. Default is 20. 
# --fastqc will run FastQC in the default mode once trimming is done. 
# TrimGalore can autodetect the adaptor if this parameter is not supplied. However, --illumina will detect any illumina adaptors.
# To control stringency of the adaptor removal specify the minimum number of required overlap with adaptor sequence. 1 is extremely strigent.
# --paired specifies if we have single-ended or paired-ended sequences. 
# compress output files with gzip.
trim_galore --paired --output_dir $OUTDIR --quality 20 --fastqc --path_to_cutadapt cutadapt --stringency 6 --gzip $INDIR/*.fastq.gz

