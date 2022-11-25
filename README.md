# Project 3. "E. coli outbreak investigation". De novo assembly and annotation of bacterial genomes.

## Step 1. Prepare the workspace and data. :shipit:

Create a directory for practice assignments moved to it. :shipit:
```ruby
mkdir Progect_3
cd Progect_3
```

### Download the raw data :shipit:
SRR292678 - paired end, insert size 470 bp 

SRR292862 – mate pair, insert size 2 kb, 

SRR292770 – mate pair, insert size 6 kb,


```ruby
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz
gzip -dk SRR292678sub_S1_L001_R1_001.fastq.gz SRR292678sub_S1_L001_R2_001.fastq.gz SRR292770_S1_L001_R1_001.fastq.gz SRR292770_S1_L001_R2_001.fastq.gz SRR292862_S2_L001_R1_001.fastq.gz SRR292862_S2_L001_R2_001.fastq.gz
```
The number of reads:

| Command                                     | Count of rows | Count of reads |
|---------------------------------------------|---------------|----------------|
| wc -l SRR292678sub_S1_L001_R1_001.fastq     | 21997384      | 5 499 346      |
| wc -l SRR292678sub_S1_L001_R2_001.fastq     | 21997384      | 5 499 346      |
| wc -l SRR292770_S1_L001_R1_001.fastq        | 20408164      | 5 102 041      |
| wc -l SRR292770_S1_L001_R2_001.fastq        | 20408164      | 5 102 041      |
| wc -l SRR292862_S2_L001_R1_001.fastq        | 20408164      | 5 102 041      |
| wc -l SRR292862_S2_L001_R2_001.fastq        | 20408164      | 5 102 041      |

### Run FastQC on all 6 fastq files. :shipit:

```ruby
CONDA_SUBDIR=osx-64 conda create -n myenv_fastqc python=3.9
conda activate myenv_fastqc
fastqc -o . SRR292678sub_S1_L001_R1_001.fastq SRR292678sub_S1_L001_R2_001.fastq
fastqc -o . SRR292770_S1_L001_R1_001.fastq SRR292770_S1_L001_R2_001.fastq
fastqc -o . SRR292862_S2_L001_R1_001.fastq SRR292862_S2_L001_R2_001.fastq
```
All data is of good quality.

## Step 1. K-mer profile and genome size estimation. :shipit:
```ruby
apt install jellyfish
jellyfish count -m 31 -o jel_outp_31_paired -C -s 50000000 SRR292678sub_S1_L001_R1_001.fastq SRR292678sub_S1_L001_R2_001.fastq
jellyfish histo -o 31_out.histo el_outp_31_paired
```
Data visualization with R:
```ruby
plot(read.table("31_out.histo")[4:600,], type="l")
```
![afd5b79f-0785-4376-bfea-aca70ea7e526](https://user-images.githubusercontent.com/109213422/203651537-adf3fec1-4c67-4a8a-946d-7854a86d0fde.png)

Knowing the number of reads and average length, estimate the total number of bases in all the reads.

**N = (M*L)/(L-K+1)**

**Genome_size = T/N**

**(N: Depth of coverage, M: Kmer peak, K: Kmer-size, L: avg read length T: Total bases)**
> - M: Kmer peak 62 
> - K: Kmer-size 31, 
> - L: avg read length 90
> - T: Total bases -- 5499346 * 90 = 494941140
> - N = (M*L)/(L-K+1) = (62 * 90) / (90 - 31 + 1) = 93
> - Genome_size = T/N = 494941140 / 93 = 5321947.7419354


## Step 3. Assembling E. coli X genome from paired reads. :shipit:

Create a directory for the processedData in Project_3.
```ruby
mkdir processedData
``` 

Run [spades](https://cab.spbu.ru/software/spades/) on short reeds and then on all reeds and compare the results with **quast** 
```ruby
spades --pe1-1 SRR292678sub_S1_L001_R1_001.fastq --pe1-2 SRR292678sub_S1_L001_R2_001.fastq -o . /processedData/spadesOutput1
spades --pe1-1 SRR292678sub_S1_L001_R1_001.fastq --pe1-2 SRR292678sub_S1_L001_L2_001.fastq --mp1-1 SRR292770_S1_L001_R1_001. fastq --mp1-2 SRR292770_S1_L001_R2_001.fastq --mp2-1 SRR292862_S2_L001_R1_001.fastq --mp2-2 SRR292862_S2_L001_R2_001.fastq -o ../processedData/spadesOutput2
```
Install [quast](https://cab.spbu.ru/software/quast/) 

```ruby
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz 
tar -xzf quast-5.2.0.tar.gz 
cd quast-5.2.0
./setup.py install
```

Run quast on the results of the previous step 
```ruby
python quast.py -o ~/IB/practice/Progect_3/processedData/quastOutput1 -s ~/IB/practice/Progect_3/processedData/spadesOutput1/scaffolds.fasta ~/IB/practice/Progect_3/processedData/spadesOutput1/contigs.fasta
python quast.py -o ~/IB/practice/Progect_3/processedData/quastOutput2 -s ~/IB/practice/Progect_3/processedData/spadesOutput2/scaffolds.fasta ~/IB/practice/Progect_3/processedData/spadesOutput2/contigs.fasta
```
![photo_2022-11-24_01-41-26](https://user-images.githubusercontent.com/109213422/203658175-84692632-d486-4408-89a8-a512aae61ebc.jpg)
![photo_2022-11-24_01-41-22](https://user-images.githubusercontent.com/109213422/203658182-15defa5d-a3f5-4799-9a65-310041ba37f8.jpg)
![photo_2022-11-24_01-41-23](https://user-images.githubusercontent.com/109213422/203658189-a1671e8f-59d7-4814-acc2-f3baa033b06d.jpg)
![photo_2022-11-24_01-41-23 (2)](https://user-images.githubusercontent.com/109213422/203658199-65ee4a67-c527-4401-a71c-b36ffb62f4f6.jpg)


## Step 4. Genome Annotation. :shipit:

Install [prokka](https://github.com/tseemann/prokka.git):
```ruby
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
apt-get update --fix-missing && apt-get install libmysqlclient21
sudo cpan Bio::Perl
git clone https://github.com/tseemann/prokka.git $HOME/prokka
apt install prokka
prokka --setupdb
```
Create new enviroment for work:
```ruby
conda create -n prokka_env -c conda-forge -c bioconda prokka
conda activate prokka_env
```

To annotate and predict genes:
```ruby
prokka --outdir ./prokka_output_my --centre X --compliant ~/IB/practice/Progect_3/processedData/spadesOutput2/scaffolds.fasta
```

## Step 5. Finding the closest relative of E. coli X.
Install Barrnap.
```ruby
git clone https://github.com/tseemann/barrnap.git
cd barrnap/bin
```
Find the 16S rRNA  in the results.
```ruby
./barrnap --quiet ~/IB/practice/Progect_3/processedData/spadesOutput2/scaffolds.fasta
```
| NODE_184_length_223_cov_0.720238_ID_565088  | barrnap:0.9 | rRNA | 95      | 205     | 5.6e-18 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |  
|---------------------------------------------|-------------|------|---------|---------|---------|---|---|---------------------------|-----------|------|---|
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 322892  | 323002  | 2.2e-11 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 323101  | 326001  | 0       | - | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 326359  | 327896  | 0       | - | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 592162  | 592272  | 6.1e-11 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 592407  | 592517  | 2.2e-11 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 592616  | 595516  | 0       | - | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 595966  | 597503  | 0       | - | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 2500844 | 2500954 | 2.2e-11 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 2501053 | 2503953 | 0       | - | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_1_length_2815616_cov_74.3819_ID_564387 | barrnap:0.9 | rRNA | 2504403 | 2505940 | 0       | - | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 43835   | 45372   | 0       | + | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 45873   | 48630   | 0       | + | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 48983   | 49084   | 5.5e-10 | + | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 85462   | 86999   | 0       | + | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 87500   | 90257   | 0       | + | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_5_length_236041_cov_85.7779_ID_563492  | barrnap:0.9 | rRNA | 90356   | 90466   | 2.2e-11 | + | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_6_length_209194_cov_80.6603_ID_563563  | barrnap:0.9 | rRNA | 111955  | 113492  | 0       | + | . | Name=16S_rRNA;product=16S | ribosomal | RNA  |   |   |   |   |
| NODE_6_length_209194_cov_80.6603_ID_563563  | barrnap:0.9 | rRNA | 113850  | 116750  | 0       | + | . | Name=23S_rRNA;product=23S | ribosomal | RNA  |   |   |   |   |
| NODE_6_length_209194_cov_80.6603_ID_563563  | barrnap:0.9 | rRNA | 116849  | 116959  | 2.2e-11 | + | . | Name=5S_rRNA;product=5S   | ribosomal | RNA  |   |   |   |   |
| NODE_71_length_720_cov_1.1218_ID_565136   | barrnap:0.9 | rRNA | 314   | 719   | 9.8e-23 | + | . | Name=16S_rRNA;product=16S | ribosomal | RNA |
| NODE_9_length_87669_cov_86.9393_ID_563775 | barrnap:0.9 | rRNA | 84700 | 84810 | 2.2e-11 | - | . | Name=5S_rRNA;product=5S   | ribosomal | RNA |
| NODE_9_length_87669_cov_86.9393_ID_563775 | barrnap:0.9 | rRNA | 84909 | 87664 | 0       | - | . | Name=23S_rRNA;product=23S | ribosomal | RNA |

## Searching for the genome in the RefSeq database with the 16S rRNA that most closely resembles the 16S rRNA we just found. :shipit:

Open the  [NCBI BLAST](http://blast.ncbi.nlm.nih.gov) homepage  and select “Nucleotide blast”. To perform the search against complete genomes in the RefSeq database, select the “Reference Genome Database (refseq_genomes)” in the “Database” field, and Escherichia coli in the “Organism” field.

The most similar sequence is [Escherichia coli 55989 (NC_011748.1)](https://www.ncbi.nlm.nih.gov/nuccore/218693476/). Download it into rowData in fasta format called 55989.fasta.

## Step 6. Search for the cause of new pathogenic signs of the bacterium. :shipit:
Install [Mauve](http://darlinglab.org/mauve/download.html)

Open “Mauve” and select “File” → “Align with progressiveMauve...”. Press “Add sequences” and select the reference genome, then the annotated E. coli X genome and start the alignment.

### Result:
| Product               | Gene                | Location        | Close                                                                            |
|-----------------------|---------------------|-----------------|----------------------------------------------------------------------------------|
| Shiga toxin           | stxB                | 348360-3483874  | nanS_5 product - putative 9-0-acelyl-N-acetylneuraminic acid deacetylase prec    |
|                       | stxA                | 3483886-3484845 | ORF B78                                                                          |
|                       |                     |                 | Phage holin/antiholin component S                                                |
|                       |                     |                 | phage lyozyme R (EC 3.2.1.17)                                                    |
|                       |                     |                 | phage antirepressor protein                                                      |
|                       |                     |                 | Phage enddopeptidase Rz                                                          |
|                       |                     |                 | ORF  B78                                                                         |
|                       |                     |                 | tRNA-arg                                                                         |
|                       |                     |                 | tRNA-arg                                                                         |
|                       |                     |                 | tRNA-met                                                                         |
|                       |                     |                 | phage DNA adenine methylase (EC 2.1.1.72)                                        |
|                       |                     |                 | phage antitermination protein Q                                                  |
|                       |                     |                 | phage protein NinH                                                               |
|                       |                     |                 | phage recombination protein NinG                                                 |
|                       |                     |                 | Gifsy-2 prophage protein                                                         |
|                       |                     |                 | DNA primase, phage associated                                                    |
|                       |                     |                 | Putative ATP-dependent helicase                                                  |
| antibiotik resistance | Bla_1 CTX-M family  | 5195566-5196441 | tryptophan synthase (indol-salavaging) (EC 4.2.1. 122)                           |
|                       |                     |                 | ProP effector                                                                    |
|                       |                     |                 | repication initiation protein                                                    |
|                       |                     |                 | Incl1 plasmid conjugative transfer protein TraA                                  |
|                       |                     |                 | Incl1 plasmid conjugative transfer NusG-type transcription a ntitterminator TraB |
|                       |                     |                 | Incl1 plasmid conjugative transfer protein TraC                                  |
|                       |                     |                 | Incl1 plasmid conjugative transfer protein PilI                                  |
|                       |                     |                 | Incl1 plasmid conjugative transfer protein PilK, M, N, O, P, R, S                |
|                       |                     |                 | Toxin coregulated pilus biosynthesis protein E                                   |
|                       |                     |                 | Transposon Tn3 resolvase                                                         |
|                       |                     |                 | class A beta-lactamase (EC 3.5.2.6) => TEM family                                |
|                       |                     |                 | errr-prone repair protein UmuD                                                   |
|                       |                     |                 | adenine-specific methyltransferase                                               |
|                       |                     |                 | antirestriction protein KlcA                                                     |
|                       |                     |                 | UPF0401 protein YkfF                                                             |
|                       |                     |                 | chromosome-pairtitioning protein spoOJ                                           |
|                       |                     |                 | CcdB toxin protein                                                               |
|                       | bla_2 TEM family    | 5199263-5200123 | transposon Tn3 resolvase                                                         |
|                       |                     |                 | class a beta-lactamase => ctx-m family, extende d-spectrum                       |
|                       |                     |                 | tryptophan synthase                                                              |
|                       |                     |                 | error-prone, lesion bypass DNA polymerase V (UmuC)                               |
|                       |                     |                 | adenine-specific methyltransferase                                               |
|                       |                     |                 | plasmid SOS inhibition protein PsiB PsiA                                         |
|                       |                     |                 | Eaa protein                                                                      |
|                       |                     |                 | 


## ResFinder result:
| Reference genome |                       |                              |                         |                    |
|------------------|-----------------------|------------------------------|-------------------------|--------------------|
|                  | Antimicrobial         | Class                        | WGS-predicted phenotype | Genetic background |
|                  | tetracycline          | tetracycline                 | Resistant               | tet(B)             |
|                  | doxycycline           | tetracycline                 | Resistant               | tet(B)             |
|                  | minocycline           | tetracycline                 | Resistant               | tet(B)             |
| Assembled genome |                       |                              |                         |                    |
|                  | cefepime              | beta-lactam                  | Resistant               | blaCTX-M-15        |
|                  | ampicillin            | beta-lactam                  | Resistant               | blaTEM-1B          |
|                  | cefotaxime            | beta-lactam                  | Resistant               | blaCTX-M-15        |
|                  | sulfamethoxazole      | folate pathway antagonist    | Resistant               | sul1               |
|                  | trimethoprim          | folate pathway antagonist    | Resistant               | dfrA7              |
|                  | tetracycline          | tetracycline                 | Resistant               | tet(A)             |
|                  | ceftazidime           | beta-lactam                  | Resistant               | blaCTX-M-15        |
|                  | streptomycin          | aminoglycoside               | Resistant               | aph(6)-Id          |
|                  | cephalothin           | beta-lactam                  | Resistant               | blaTEM-1B          |
|                  | piperacillin          | beta-lactam                  | Resistant               | blaTEM-1B          |
|                  | amoxicillin           | beta-lactam                  | Resistant               | blaTEM-1B          |
|                  | ethidium bromide      | quaternary ammonium compound | Resistant               | qacE               |
|                  | doxycycline           | tetracycline                 | Resistant               | tet(A)             |
|                  | benzylkonium chloride | quaternary ammonium compound | Resistant               | qacE               |
|                  |                       |                              |                         |                    |
|                  | ticarcillin           | beta-lactam                  | Resistant               | blaTEM-1B          |
|                  | chlorhexidine         | quaternary ammonium compound | Resistant               | qacE               |
