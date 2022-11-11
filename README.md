# Progect 2

### Create a directory and environment for practice assignments moved to it. 

```ruby
CONDA_SUBDIR=osx-64 conda create -n myenv_2 python=3.9
conda activate myenv_2
mkdir practice/Progect_2
cd practice/Progect_2
```


### Let's install the utilities we need for our work.
```ruby
sudo apt update
sudo apt install bwa samtools bedtools snakemake
```
### Download the study data and reference for the reference sequence of the parental E.coli strain
```ruby
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz
```
[Reference](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta): `select Send to -> File -> FASTA format`. It had been saved as `sequence.fasta`

Unpack the fastq files 
```ruby
gzip -dk SRR1705851.fastq.gz
```
### Run fastqc
```ruby
fastqc -o . SRR1705851.fastq
```
All of our readings were of satisfactory quality, so we don't need to do the trimming.

```ruby
conda deactivate
```
###Work with fasta file

Create a genomic index by adding -p flag to the index name.    
```ruby
bwa index sequence.fasta -p index.bwa_index
```
Aligh the reads to the resulting genome index and convert a **sam** file to a **bam** file
```ruby
bwa mem -t 4 index.bwa_index SRR1705851.fastq | samtools view -S -b -
```
Counts the number of alignments for each FLAG type.
```ruby
samtools flagstat bwa_output.bam
```
**Output:**
> - 361349 + 0 in total (QC-passed reads + QC-failed reads)
> - 0 + 0 secondary
> - 3084 + 0 supplementary
> - 0 + 0 duplicates
> - 361116 + 0 mapped (99.94% : N/A)
> - 0 + 0 paired in sequencing
> - 0 + 0 read1
> - 0 + 0 read2
> - 0 + 0 properly paired (N/A : N/A)
> - 0 + 0 with itself and mate mapped
> - 0 + 0 singletons (N/A : N/A)
> - 0 + 0 with mate mapped to a different chr
> - 0 + 0 with mate mapped to a different chr (mapQ>=5)

Sort and index bam file
```ruby
samtools sort bwa_output.bam -o | samtools index -
```
Index of alignment files and obtaining mpileup files, the depth limit in mpileup was removed (parameter -d 0).
```ruby
samtools mpileup -d 0 -f sequence.fasta bwa_output_sorted.bam > my_2.mpileup
```
*We take all reads, since our coverage there is large enough, entering the -d 0 parameter allows us to find a rare mutation, and usually once the number of reads reaches 8000, they stop being counted and scanned.*

Find rare variants in sequencing data
```ruby
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp my_2.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_2.vcf
```
**Output:**
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.001
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from my_2.mpileup
> - 1665 bases in pileup file
> - 23 variant positions (21 SNP, 2 indel)
> - 0 were failed by the strand-filter
> - 21 variant positions reported (21 SNP, 0 indel)

Look results
```ruby
cat VarScan_results_001.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > rommate.txt
```
**Output:**

|Name | Position | Original nucleotide| Replacement | FREQ |
|-----------|----|---|---|--------------------------------------|
|KF848938.1 | 72 | A | G | 1/1:255:16832:16794:6:16787:99.96%:0E0:35:36:4:2:10898:5889|
|KF848938.1 | 117 | C | T | 1/1:255:20768:20663:36:20625:99.82%:0E0:35:37:27:9:13462:7163
|KF848938.1 | 254 | A | G | 0/1:20:35781:35626:35562:59:0.17%:8.5683E-3:36:36:23919:11643:37:22
|KF848938.1 | 276 | A | G | 0/1:24:37022:36965:36892:64:0.17%:3.3004E-3:37:35:22579:14313:30:34
|KF848938.1 | 307 | C | T | 0/1:255:37506:37386:37029:351:0.94%:6.9068E-66:36:35:22400:14629:184:167
|KF848938.1 | 340 | T | C | 0/1:23:37973:37793:37723:64:0.17%:4.6441E-3:37:36:23413:14310:40:24

###Download fastq data for the three controls (from sequencing of isogenic reference samples) from SRA FTP:
```ruby
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz # SRR1705858
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz # SRR1705859
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz # SRR1705860
```

### Look at the number of lines and count the number of reads in each file
```ruby
zcat SRR1705858.fastq.gz | wc -l
```
`1026344 / 4  = 256 586`

```ruby
zcat SRR1705859.fastq.gz | wc -l
```
`933308 / 4 = 233 327`

```ruby
zcat SRR1705860.fastq.gz | wc -l
```
`999856 / 4 = 249 964`

###Do the same preparing:
```ruby
bwa mem -t 4 index.bwa_index SRR1705858.fastq > bwa_output_SRR1705858.sam
bwa mem -t 4 index.bwa_index SRR1705859.fastq > bwa_output_SRR1705859.sam
bwa mem -t 4 index.bwa_index SRR1705860.fastq > bwa_output_SRR1705860.sam
 
samtools view -S -b bwa_output_SRR1705858.sam > bwa_output_SRR1705858.bam
samtools view -S -b bwa_output_SRR1705859.sam > bwa_output_SRR1705859.bam
samtools view -S -b bwa_output_SRR1705860.sam > bwa_output_SRR1705860.bam
 
samtools sort bwa_output_SRR1705858.bam -o bwa_output_sorted_SRR1705858.bam
samtools sort bwa_output_SRR1705859.bam -o bwa_output_sorted_SRR1705859.bam
samtools sort bwa_output_SRR1705860.bam -o bwa_output_sorted_SRR1705860.bam

samtools index bwa_output_sorted_SRR1705858.bam
samtools index bwa_output_sorted_SRR1705859.bam
samtools index bwa_output_sorted_SRR1705860.bam

samtools mpileup -d 0 -f sequence.fasta bwa_output_sorted_SRR1705858.bam > my_SRR1705858.mpileup
samtools mpileup -d 0 -f sequence.fasta bwa_output_sorted_SRR1705859.bam > my_SRR1705859.mpileup
samtools mpileup -d 0 -f sequence.fasta bwa_output_sorted_SRR1705860.bam > my_SRR1705860.mpileup
```
```ruby
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp my_SRR1705858.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_SRR1705858.vcf
```
**Output:**
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.001
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from my_SRR1705858.mpileup
> - 1665 bases in pileup file
> - 58 variant positions (58 SNP, 0 indel)
> - 1 were failed by the strand-filter
> - **57 variant positions reported (57 SNP, 0 indel)**

```ruby
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp my_SRR1705859.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_SRR1705859.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.001
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from my_SRR1705859.mpileup
> - 1665 bases in pileup file
> - 54 variant positions (54 SNP, 0 indel)
> - 2 were failed by the strand-filter
> - **52 variant positions reported (52 SNP, 0 indel)**

```ruby
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp my_SRR1705860.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_SRR1705860.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.001
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from my_SRR1705860.mpileup
> - 1665 bases in pileup file
> - 61 variant positions (61 SNP, 0 indel)
> - 0 were failed by the strand-filter
> - **61 variant positions reported (61 SNP, 0 indel)**

###Check deph of reading considering all files:
```ruby
samtools coverage --reference sequence.fasta bwa_output_sorted_SRR1705858.bam bwa_output_sorted_SRR1705859.bam bwa_output_sorted_SRR1705860.bam
```
|rname | startpos | endpos | numreads | covbases | coverage | meandepth | meanbaseq | meanmapq|
|------|----------|--------|----------|----------|----------|-----------|----------|---------|
|KF848938.1  |    1    |    1665  |  740141    |    1665      |     100       |     65244.4      |     36.6     |      60|

```ruby
cat VarScan_results_SRR1705858.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > SRR1705858.txt
cat VarScan_results_SRR1705859.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > SRR1705859.txt
cat VarScan_results_SRR1705860.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > SRR1705860.txt
```
Next, for each file, we count the average mutation frequency and the average deviation. Next, in the file with the studied data, we find mutations whose frequency of occurrence deviates by more than 3 standard deviations from the mean of the three averages.
