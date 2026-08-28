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
