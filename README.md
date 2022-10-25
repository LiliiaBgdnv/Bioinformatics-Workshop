# Progect 1. “What causes antibiotic resistance?”. Alignment to reference, variant calling.

##▶The first step. Prepare the workspace.

###Create a directory for practice assignments moved to it. 

```ruby
mkdir practice/Progect_1
cd practice/Progect_1
```

###Download the raw data and annotation for the reference sequence of the parental E.coli strain

```ruby
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
```
Also we need raw Illumina equencing reads. Just download two files from [here](https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3). It is forward and reveerse reads.

##▶The second step. Let's check our data.

If the download was successful we will have 4 files:

> - `GCF_000005845.2_ASM584v2_genomic.fna.gz` — reference sequence
> - `GCF_000005845.2_ASM584v2_genomic.gff.gz` — annotation
> - `amp_res_1.fastq.gz` — forward reads
> - `amp_res_2.fastq.gz` — reverse reads

Look at them in more detail.:
```ruby
zcat amp_res_1.fastq.gz | head -4
```
Information about one reading is contained in 4 lines: 
1. A sequence identifier with information about the sequencing run and the cluster.
2. The sequence.
3. A separator, which is simply a plus (+) sign.
4. The base call quality scores. These are Phred +33 encoded characters to represent the numerical quality scores.

| Line |
|-------------|
| `@SRR1363257.37 GWZHISEQ01:153:C1W31ACXX:5:1101:14027:2198 length=101`    |
|`GGTTGCAGATTCGCAGTGTCGCTGTTCCAGCGCATCACATCTTTGATGTTCACGCCGTGGCGTACACG`     |
| `+`  |
| `@?:=:;DBFADH;CAECEE@@:FFHGAE4?C?DE<BFGEC>?>FHE4BFFIIFHIBABEECA83;>>@`     |

Check how many lines there are in each fastq file.

```ruby
zcat amp_res_1.fastq.gz | wc -l
zcat amp_res_2.fastq.gz | wc -l
```
In each file there are `1823504` lines, if we divide the number of lines by 4 (which is how many lines are allocated per read), we get `455876` reads.

##▶The third step. Inspect raw sequencing data with fastqc. Filtering the reads.

**Let's download FASTQC.**

> To do this, create a conda environment. To avoid problems with conflicting java versions of different programs in the future.
> ```ruby
> CONDA_SUBDIR=osx-64 conda create -n myenv_x86 python=3.9
> ```
>Activate the environment
>```ruby
>conda activate myenv_x86
> conda config --env --set subdir osx-64
>```
> Install fastqc
>```ruby
>apt-get install fastqc 
>```
>Check if the file is installed 
>```ruby
>fastqc -h  
>```
>
Unpack the fastq files 
```ruby
gzip -dk amp_res_1.fastq.gz
gzip -dk amp_res_2.fastq.gz
```
`flag -d indicates that we need exactly unpack, flag -k  indicates that we need to save the compressed files`

###Run fastqc
```ruby
fastqc -o . amp_res_1.fastq amp_res_2.fastq
```
For each file we get a file with the extension `.html`. Check them out.![photo_2022-10-25_15-16-04.jpg]

##▶The fourth step.  Filtering the reads.

**To filter the readings we will use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).**

To work with Trimmomatic, let's leave the environment
```ruby
conda deactivate
```
Installing and Checking
```ruby
apt-get install trimmomatic
TrimmomaticPE
```

>The parameters are summarized as follows:
- trimmomatic (**PE** if paired | **SE** if single reads) 
- **-phred33 | -phred64** (for Sanger / Illumina 1.9) 
- input.fastq output.fastq 
- **LEADING:...** (remove the quality from the beginning with a quality lower than specified) 
- **TRAILING:...** (cut off end nucleotides with quality below the specified number) 
- **MINLEN:...**(discard all reads with length less than a specified number)
- **SLIDINGWINDOW:...:...** (window width of number 1 when quality per base falls below number 2)
>

###Run Trimmomatic in paired end mode, with following parameters:
* Cut bases off the start of a read if quality below 20 (`LEADING`)
* Cut bases off the end of a read if quality below 20. (`TRAILING`)
* Trim reads using a sliding window approach, with window size 10 and average quality within the window 20. (`SLIDINGWINDOW:10:20`)
* Drop the read if it is below length 20. (`MINLEN`)
```
trimmomatic PE -phred33 amp_res_1.fastq amp_res_2.fastq amp_res_1_trimmed.fastq amp_res_1_solo.fastq amp_res_2_trimmed.fastq amp_res_2_solo.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```

After the Trimmomatic runs, we get:

`Input Read Pairs: 455876 Both Surviving: 446259 (97.89%) Forward Only Surviving: 9216 (2.02%) Reverse Only Surviving: 273 (0.06%) Dropped: 128 (0.03%)`

Let's check the number of lines for reads was saved for forward and reverse reads:
``` ruby
less amp_res_1_trimmed.fastq | wc -l
less amp_res_2_trimmed.fastq | wc -l
```
In each file there are `1785036` lines, if we divide the number of lines by 4 (which is how many lines are allocated per read), we get `446259` reads.

###Repeat the fastqc analysis
```ruby
conda activate myenv_x86
fastqc -o . amp_res_1_trimmed.fastq amp_res_2_trimmed.fastq
```
##▶The five step. Aligning sequences to reference
Unpack the archive with the reference E.coli sequence 
```ruby
gzip -dk GCF_000005845.2_ASM584v2_genomic.fna.gz
```
Download aligner
```ruby
apt-get install bwa 
```
Let's create a genomic index by adding -p flag to the index name.    
```ruby
bwa index GCF_000005845.2_ASM584v2_genomic.fna -p index.bwa_index
```
We get four files :
> - `index.bwa_index.amb`
> - `index.bwa_index.ann`
> - `index.bwa_index.bwt`
> - `index.bwa_index.pac`
> - `index.bwa_index.sa`

Aligh the reads to the resulting genome index, give the `number of threads` with the **-t** flag, `file with the indexes`, 
`file with direct reads`, `reverse reads`, and `output file` with the extension **sam** through the symbol **>**
```ruby
bwa mem -t 4 index.bwa_index trimm_part/amp_res_1_trimmed.fastq trimm_part/amp_res_2_trimmed.fastq > bwa_output.sam
```
Convert a **sam** file to a **bam** file `(the -S -b flags indicate from which format you want to convert)`
```ruby
samtools view -S -b bwa_output.sam > bwa_output.bam
```
To get some basic statistics.
```ruby
samtools flagstat bwa_output.bam
```
>- 892776 + 0 in total (QC-passed reads + QC-failed reads)
> - 0 + 0 secondary
> - 258 + 0 supplementary
> - 0 + 0 duplicates
> - 891649 + 0 mapped (99.87% : N/A)
> - 892518 + 0 paired in sequencing
> - 446259 + 0 read1
> - 446259 + 0 read2
> - 888554 + 0 properly paired (99.56% : N/A)
> - 890412 + 0 with itself and mate mapped
> - 979 + 0 singletons (0.11% : N/A)
> - 0 + 0 with mate mapped to a different chr
> - 0 + 0 with mate mapped to a different chr (mapQ>=5)

Sort bam file (the program sort sorts bam file by chromosome and by coordinate)
```ruby
- samtools sort bwa_output.bam -o bwa_output_sorted.bam
```
Index the sorted file.
```ruby
samtools index bwa_output_sorted.bam
```
##▶The six step. Variant calling
*Let's look at our data on the reference genome and how many reads in our data have mutations and in which position. To do this, we're going to need [VarScan](http://dkoboldt.github.io/varscan/) (variant scanner).*

###Make a sorted and indexed bam file
```ruby
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna bwa_output_sorted.bam > my.mpileup
```
###Installing and checking VarScan
```ruby
conda install -c bioconda varscan
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp -h
```
###Let's run the program with a 50% threshold 
The option --variants flag tells VarScan to output only positions that
above our threshold, and the option --output-vcf 1 specifies that we want
receive data in another kind of data format, called vcf (variant call format)
```ruby
java -jar ~/miniconda3/share/varscan-2.4.4-1/VarScan.jar mpileup2snp my.mpileup --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.5
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from GCF_000005845.2_ASM584v2_genomic.mpileup
> - 4641343 bases in pileup file
> - 9 variant positions (6 SNP, 3 indel)
> - 1 were failed by the strand-filter
> - 5 variant positions reported (5 SNP, 0 indel)

##▶The seven step. Automatic SNP annotation
To avoid making unnecessary mistakes in the annotation, we will use the automatic tool for snp annotation. For this project we can use [SnpEff](http://pcingola.github.io/SnpEff/) (short for “SNP effect”).

### Installing SnpEff
```ruby
conda install -c bioconda snpEff
```
### Download reference sequence 

```ruby
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz
```
### Create database

> Create folder for the database
```ruby
mkdir -p data/k12
```
Unzip the reference sequence and move it to the folder with the name change
```ruby
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
cp GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk
```
> 
Create database
```ruby
snpEff build -genbank -v k12
```
### Annotation
```ruby
snpEff ann k12 VarScan_results.vcf >  VarScan_results_annotated.vcf
```
##▶The eight step. Visualization.😎
Now it's time to see  with our eyes the data we processed and the mutations we found. Do it in the [IGV](https://software.broadinstitute.org/software/igv/).
You may use online or localy on your computer.

> Select “Genomes”, ”Load Genome from File” and select **GCF_000005845.2_ASM584v2_genomic.fna**.
Then select “File”, “Open from file” and select **bwa_output_sorted.bam,  VarScan_results.vcf,  VarScan_results_annotated.vcf,  GCF_000005845.2_ASM584v2_genomic.gff**.
