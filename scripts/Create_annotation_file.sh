INDIR="Genome"
OUTDIR="Genome/Stringtie"
REF="vu-2k.fna" # Adjust depending on the genome needed

ml Miniforge3

source activate liftoff 

liftoff -g $INDIR/GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.gtf \
	-o $INDIR/liftoff_vu2k.gtf \
	-p 8 -a 0.9 -s 0.9 \
	$INDIR/vu-2k.fasta $INDIR/GCF_900497805.2_bare-nosed_wombat_genome_assembly_genomic.fna
