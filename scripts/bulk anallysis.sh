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

