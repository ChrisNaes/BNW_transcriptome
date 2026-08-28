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

