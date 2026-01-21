#!/bin/bash
#SBATCH --job-name=FeatureCounts				# Job name
#SBATCH --partition=batch						# Partition (queue) name
#SBATCH --ntasks=100							# Run on a single CPU
#SBATCH --mem=20gb								# Job memory request
#SBATCH --time=100:00:100						# Time limit hrs:min:sec
#SBATCH --output=/home/cn68176/Log/%x.%j.out    # Standard output log
#SBATCH --error=/home/cn68176/Log/%x.%j.err		# Standard error log
#SBATCH --mail-type=END,FAIL          			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cn68176@uga.edu   			# Where to send mail

INDIR="/scratch/cn68176/RNAseq_BNW/Alignment"
OUTDIR="/scratch/cn68176/RNAseq_BNW/GeneCounts"
REFDIR="/scratch/cn68176/RNAseq_BNW/Genome" # Adjust depending on the genome needed
REF="vu-2k.fna" # Adjust depending on the genome needed
REFANN="/scratch/cn68176/RNAseq_BNW/Genome/Stringtie_wRef/Thyroid_fixed_gene_id.gtf"

ml SAMtools/1.18-GCC-12.3.0
ml Subread/2.0.6-GCC-12.3.0

## We first use samtools to create sorted and indexed BAM files for our samples
#for FILE in $INDIR/*.bam; do
	#Extract the base filename 
#	base_name=$(basename "$FILE" .bam)
	# Convert SAM to BAM files
#	samtools view -bS "$INDIR/${base_name}.sam" > "$INDIR/${base_name}.bam"

	# Sort BAM files
#	samtools sort -@ 8 "$INDIR/${base_name}.bam" -o "$INDIR/${base_name}_sorted.bam"

	# Index the sorted BAM files
#	samtools index "$INDIR/${base_name}_sorted.bam"

	# Print the information for the sample that completes then remove the temporary files
#	echo ${base_name}.bam
	#rm $INDIR/${base_name}.bam
	#rm $INDIR/${base_name}.sam

#done

#######################################
#samtools merge -o Ant.bam Ant_S8_L001.bam Ant_S8_L002.bam Ant_S8_L003.bam Ant_S8_L004.bam
#samtools merge -o Bon.bam Bon_S1_L001.bam Bon_S1_L002.bam Bon_S1_L003.bam Bon_S1_L004.bam
#samtools merge -o Both.bam Both_S5_L001.bam Both_S5_L002.bam Both_S5_L003.bam Both_S5_L004.bam
#samtools merge -o Bri.bam Bri_S10_L001.bam Bri_S10_L002.bam Bri_S10_L003.bam Bri_S10_L004.bam
#samtools merge -o Cle.bam Cle_S7_L001.bam Cle_S7_L002.bam Cle_S7_L003.bam Cle_S7_L004.bam
#samtools merge -o Col.bam Col_S6_L001.bam Col_S6_L002.bam Col_S6_L003.bam Col_S6_L004.bam
#samtools merge -o DW01.bam DW01_S9_L001.bam DW01_S9_L002.bam DW01_S9_L003.bam DW01_S9_L004.bam
#samtools merge -o DW02.bam DW02_S2_L001.bam DW02_S2_L002.bam DW02_S2_L003.bam DW02_S2_L004.bam
#samtools merge -o Nub.bam Nub_S3_L001.bam Nub_S3_L002.bam Nub_S3_L003.bam Nub_S3_L004.bam
#samtools merge -o Vet.bam Vet_S4_L001.bam Vet_S4_L002.bam Vet_S4_L003.bam Vet_S4_L004.bam

# Then we use featureCounts to quantify expressed genes and transcripts
	# -T is number of threads 
	# -g is 
	# -a is your annotation (GTF) file
	# -o is your output file, in this case I'm saving it as a txt file.
	# last command is providing your sorted and indexed BAM file.
featureCounts -T 32 -t exon -g gene_id \
	-a $REFANN \
	-o ${OUTDIR}/Thyroid.txt \
	$INDIR/*_sorted.bam




#featureCounts -T 4 -a $REFDIR/test_vour.gtf -o test_counts.txt -t exon -g gene_id Ant_sorted.bam


