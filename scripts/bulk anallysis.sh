fastqc \
  ~/data/rnaseq/raw_fastq/*.fastq \
  ~/data/rnaseq/raw_fastq/*.fastq.gz \
  --outdir ~/data/rnaseq/qc


# Get Gene Fasta (primary assembly)
mkdir -p ~/ref/ensembl_114/{00_raw,01_fasta,02_gtf}
cd ~/ref/ensembl_114/00_raw

wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  > ../01_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Gene annotation GTF
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

gunzip -c Homo_sapiens.GRCh38.114.gtf.gz \
  > ../02_gtf/Homo_sapiens.GRCh38.114.gtf


# ----- trim 14 bp----

IN=~/data/rnaseq/raw_fastq
TRIMIN=~/data/rnaseq/trimmed_fastq
mkdir -p "$TRIMIN"

for fq in "$IN"/*.fastq.gz "$IN"/*.fq.gz; do
  [ -e "$fq" ] || continue

  base=$(basename "$fq")
  sample=${base%.fastq.gz}
  sample=${sample%.fq.gz}

  cutadapt -u 14 -o "$TRIMIN/${sample}.trim14.fastq.gz" "$fq"
done


#Create index
GBASE=~/ref/ensembl_114
FA=$GBASE/01_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=$GBASE/02_gtf/Homo_sapiens.GRCh38.114.gtf
GDIR=$GBASE/10_star_index

mkdir -p "$GDIR"

~/tools/STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir "$GDIR" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 93


#Save index with specific overhand
# Create index (sjdbOverhang = 79)
GBASE=~/ref/ensembl_114
FA=$GBASE/01_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=$GBASE/02_gtf/Homo_sapiens.GRCh38.114.gtf
GDIR=$GBASE/10_star_index_79

mkdir -p "$GDIR"

~/tools/STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir "$GDIR" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 79



## Align trimmed

GBASE=~/ref/ensembl_114
GDIR=$GBASE/10_star_index_79

TRIMIN=~/data/rnaseq/trimmed_fastq
OUTBASE=~/data/rnaseq/aligned_trim14
mkdir -p "$OUTBASE"

for fq in "$TRIMIN"/*.trim14.fastq.gz; do
  [ -e "$fq" ] || continue

  base=$(basename "$fq")
  sample=${base%.trim14.fastq.gz}

  outdir="$OUTBASE/$sample"
  mkdir -p "$outdir"

  ~/tools/STAR --runThreadN 8 \
    --genomeDir "$GDIR" \
    --readFilesIn "$fq" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$outdir/" \
    --outSAMtype BAM SortedByCoordinate
done




## Count reads
GTF=~/ref/ensembl_114/02_gtf/Homo_sapiens.GRCh38.114.gtf
ALIGNED=~/data/rnaseq/aligned_trim14
QCBASE=~/data/rnaseq/qc/featureCounts_trim14
THREADS=20

mkdir -p "$QCBASE"

for BAM in "$ALIGNED"/*/Aligned.sortedByCoord.out.bam; do
  [ -e "$BAM" ] || continue

  SAMPLE=$(basename "$(dirname "$BAM")")
  OUTDIR="$QCBASE/$SAMPLE"
  mkdir -p "$OUTDIR"

  featureCounts -T "$THREADS" \
    -t exon \
    -g gene_id \
    -s 1 \
    -a "$GTF" \
    -o "$OUTDIR/${SAMPLE}_featureCounts_exon.txt" \
    "$BAM"
done



