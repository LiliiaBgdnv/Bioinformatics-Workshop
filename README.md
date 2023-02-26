# Progect 6. RNA-seq. Baking Bread.

New bioinformatics skills covered: splice-junction aware alignment,
guided transcript assembly, differential expression analysis

## Step 1. Download data.
```ruby
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941816/SRR941816.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941817/SRR941817.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941818/SRR941818.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941819/SRR941819.fastq.gz
```

As a reference genome we will use Saccharomyces cerevisiae:
```ruby 
wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
```

## Step 2. Analysis.

Install [hisat2](https://daehwankimlab.github.io/hisat2/)
```ruby 
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
make
```

build genome index:
```ruby
hisat2-build GCF_000146045.2_R64_genomic.fna GCF_000146045.2_R64_genomic

./hisat2/hisat2 -p 2 -x GCF_000146045.2_R64_genomic -U SRR941816.fastq | samtools sort > out16.bam
./hisat2/hisat2 -p 2 -x GCF_000146045.2_R64_genomic -U SRR941817.fastq | samtools sort > out17.bam
./hisat2/hisat2 -p 2 -x GCF_000146045.2_R64_genomic -U SRR941818.fastq | samtools sort > out18.bam
./hisat2/hisat2 -p 2 -x GCF_000146045.2_R64_genomic -U SRR941819.fastq | samtools sort > out19.bam
```
Quantifying with featureCounts, first of all convert the GFF file to GTF format:
```ruby
conda install gffread 
gffread GCF_000146045.2_R64_genomic.gff -T -o GCF_000146045.2_R64_genomic.gft
featureCounts -g gene_id -a GCF_000146045.2_R64_genomic.gft -o result out16.bam out17.bam out18.bam out19.baml
cat result | cut -f 1,7-10 > simple_counts.txt
```
## Step 3. Find differentially expressed genes with Deseq2.
Change in R script:
```ruby
#add:
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
library(locfit)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
```

calculate metrics:
```ruby
cat simple_counts.txt | R -f deseq2.r
cat norm-matrix-deseq2.txt | R -f draw-heatmap.r
head -n 50 result.txt | cut -f 1 | cut -d "-" -f 2 > genes.txt
```

## Step 4.  [Result Interpretation](http://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl).

For our top 50 differentially expressed genes:
- in step 1 press “Choose file” and upload genes.txt 
- in step 2, select “Yeast GO-Slim: Process”
- in step 3, make sure “SELECT ALL TERMS” is highlighted. Press “Search”
- Try to interpret these [results](https://gotermfinder.yeastgenome.org/mapper_genes_1456146_slimTerms.html)

