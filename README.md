# Project 4. Tardigrades: from genestealers to space marines :pig2:

*New bioinformatics skills covered: gene prediction, functional annotation, protein localization*

Tardigrades are tiny creatures that can survive in some of the worst environmental circumstances. They are also known as water bears, chubby wudgies, or moss piglets. The deep sea to the Himalayas, these eight-legged, aquatic animals may be found all over the planet. They are referred to as "extremophiles" since they can endure radiation, complete dehydration, pressure, and cold.

The ability to enter a dehydrated condition (and remain there for up to five years) and the fact that this state offers fewer reactants for ionizing radiation were used by researchers to explain radiation resistance. But further research revealed that, compared to other animals, they nevertheless have a high level of resistance to shortwave UV light when hydrated. 

It implies that they are capable of effectively repairing DNA damage. Up until 2016, when we were able to sequence the genome of Ramazzottius varieornatus, one of the Tardigrade species that is most stress-tolerant. We will try to examine this genome and uncover this secret in this research.

## Step 1. Obtaining data. Genome sequence. :pig2:
```ruby
mkdir Progect_4
cd Progect_4
```
Download and unzip assembled genome here: 
```ruby
wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz
gzip -dk GCA_001949185.1_Rvar_4.0_genomic.fna.gz
```

## Step 2. Functional annotation. :pig2:
Script:
```ruby
wget http://augustus.gobics.de/binaries/scripts/getAnnoFasta.pl
```
Extract protein sequences (fasta) from the prediction output.
```ruby
perl getAnnoFasta.pl augustus.whole.gff
```
Count the number of obtained proteins:
```ruby
grep '>' augustus.whole.aa | wc -l
```
Result: **16435**

## Step 3. Physical localization. :pig2:
Do a local alignment-based search: create a local database from your protein fasta file and look it up using your peptide sequence file as a query.

### Perform the search:

**Way a:**
```ruby
conda install -c bioconda blast
```
Create a database: 
```ruby
makeblastdb -in augustus.whole.aa -dbtype prot -out tardigrades_db
```
Search:
```ruby
blastp -db tardigrades_db -query peptides.fa -outfmt 6 -out protein_and_pept
ide.tsv
```
**Way b:**

Extract proteins of interest from the initial file.
```ruby
conda install -c bioconda diamond
```
Create a database: 
```ruby
mkdir diamond
diamond makedb --in augustus.whole.aa --db diamond/tardigrades_db
```
Search:
```ruby
diamond blastp -d diamond/tardigrades_db -q peptides.fa -f 6 -o protein_and_peptide.tsv
```
Result using python:
```ruby
names = ['peptide_id', 'protein_id', 'match_percentage', 'length', 
         'n_mismatches', 'n_gap_opening', 'start_pos_of_peptide',        'end_pos_of_peptide', 'start_pos_of_protein', 'end_pos_of_protein', 'evalue', 'bitscore']
protein_and_peptide = pd.read_csv('protein_and_peptide.tsv', 
                                            sep='\t', names=names)
protein_and_peptide
```
| peptide_id | protein_id | match_percentage | length  | n_mismatches | n_gap_opening | start_pos_of_peptide | end_pos_of_peptide | start_pos_of_protein | end_pos_of_protein | evalue | bitscore |
|------------|------------|------------------|---------|--------------|---------------|----------------------|--------------------|----------------------|--------------------|--------|----------|
| 0          | 1          | g5641.t1         | 100.000 | 9            | 0             | 0                    | 1                  | 9                    | 25                 | 33     | 2.000    | 21.6 |
| 1          | 1          | g15153.t1        | 100.000 | 9            | 0             | 0                    | 1                  | 9                    | 25                 | 33     | 2.100    | 21.6 |
| 2          | 2          | g5641.t1         | 100.000 | 9            | 0             | 0                    | 1                  | 9                    | 36                 | 44     | 1.500    | 21.9 |
| 3          | 2          | g12562.t1        | 88.889  | 9            | 1             | 0                    | 1                  | 9                    | 36                 | 44     | 3.300    | 20.8 |
| 4          | 2          | g5616.t1         | 88.889  | 9            | 1             | 0                    | 1                  | 9                    | 36                 | 44     | 6.700    | 20.0 |
| ...        | ...        | ...              | ...     | ...          | ...           | ...                  | ...                | ...                  | ...                | ...    | ...      | ...  |
| 113        | 40         | g5467.t1         | 88.889  | 9            | 1             | 0                    | 1                  | 9                    | 113                | 121    | 5.300    | 20.4 |
| 114        | 41         | g3428.t1         | 100.000 | 11           | 0             | 0                    | 1                  | 11                   | 160                | 170    | 0.055    | 26.2 |
| 115        | 41         | g11513.t1        | 87.500  | 8            | 1             | 0                    | 2                  | 9                    | 872                | 879    | 8.300    | 20.0 |
| 116        | 42         | g3679.t1         | 100.000 | 12           | 0             | 0                    | 1                  | 12                   | 21                 | 32     | 0.046    | 26.6 |
| 117        | 43         | g5237.t1         | 100.000 | 13           | 0             | 0                    | 1                  | 13                   | 91                 | 103    | 0.004    | 29.6 |

Extract proteins of interest from the initial file:
```ruby
xargs samtools faidx augustus.whole.aa < proteins_of_interest_no_duplications.txt > proteins_of_interest.faa
```
## [Result.](output_protein_nucl.txt)

## Step 4. Localization prediction. :pig2:

[WoLF PSORT](https://wolfpsort.hgc.jp/)

| protein_id | localization | e-value |
|------------|--------------|-------|
| g16368.t1  | nucl         | 20.5  |
| g10514.t1  | nucl         | 19    |
| g14472.t1  | nucl         | 28    |
| g16318.t1  | nucl         | 20.5  |
| g8312.t1   | nucl         | 15.5  |
| g15484.t1  | nucl         | 17.5  |
| g7861.t1   | nucl         | 16    |
| g11806.t1  | nucl         | 18    |
| g11960.t1  | nucl         | 32    |
| g8100.t1   | nucl         | 16.5  |
| g5927.t1   | nucl         | 30.5  |
| g10513.t1  | nucl         | 20    |

[TargetP 1.1](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)

| ID        | Prediction | OTHER    | SP       | mTP      | CS Position |
|-----------|------------|----------|----------|----------|-------------|
| g16368.t1 | OTHER      | 0.996693 | 0.003307 | 0.000000 | NaN         |
| g10514.t1 | OTHER      | 0.999543 | 0.000349 | 0.000107 | NaN         |
| g4106.t1  | OTHER      | 0.729658 | 0.266917 | 0.003425 | NaN         |
| g4970.t1  | OTHER      | 0.999996 | 0.000003 | 0.000001 | NaN         |
| g5510.t1  | OTHER      | 0.999108 | 0.000016 | 0.000876 | NaN         |
| g14472.t1 | OTHER      | 0.999999 | 0.000001 | 0.000000 | NaN         |
| g16318.t1 | OTHER      | 0.997047 | 0.002953 | 0.000000 | NaN         |
| g5237.t1  | OTHER      | 0.999545 | 0.000345 | 0.000111 | NaN         |
| g8312.t1  | OTHER      | 0.999930 | 0.000065 | 0.000004 | NaN         |
| g12510.t1 | OTHER      | 0.999738 | 0.000099 | 0.000163 | NaN         |
| g15484.t1 | OTHER      | 0.999980 | 0.000010 | 0.000010 | NaN         |
| g2203.t1  | OTHER      | 0.999869 | 0.000031 | 0.000100 | NaN         |
| g3428.t1  | OTHER      | 0.999903 | 0.000033 | 0.000064 | NaN         |
| g5443.t1  | OTHER      | 0.952853 | 0.043784 | 0.003363 | NaN         |
| g7861.t1  | OTHER      | 0.999975 | 0.000004 | 0.000022 | NaN         |
| g11806.t1 | OTHER      | 0.998977 | 0.000887 | 0.000136 | NaN         |
| g11513.t1 | OTHER      | 0.999434 | 0.000401 | 0.000164 | NaN         |
| g11960.t1 | OTHER      | 0.999996 | 0.000002 | 0.000002 | NaN         |
| g8100.t1  | OTHER      | 0.999955 | 0.000024 | 0.000021 | NaN         |
| g5927.t1  | OTHER      | 0.999995 | 0.000001 | 0.000004 | NaN         |
| g10513.t1 | OTHER      | 0.999999 | 0.000001 | 0.000000 | NaN         |

## Step 5. [BLAST](https://blast.ncbi.nlm.nih.gov) search. :pig2:

| Protein_id | DB | Accession   | Description                                                                                         | Organism                  | Length | Score(Bits) | Identities(%) | Positives(%) | E()      |
|------------|----|-------------|-----------------------------------------------------------------------------------------------------|---------------------------|--------|-------------|---------------|--------------|----------|
| g16318.t1  | TR | A2VD00.1    | Eukaryotic translation initiation factor 3 subunit A                                                | Xenopus laevis            | 1424   | 36.1        | 36.11         | 57.4         | 4,00E-08 |
| g14472.t1  | TR | P0DOW4.1    | Damage suppressor protein [Ramazzottius varieornatus]                                               | Ramazzottius varieornatus | 445    | 814         | 100           | 100          | 0        |
| g5927.t1   | TR | Q17427.1    | Glucosamine 6-phosphate N-acetyltransferase                                                         | Caenorhabditis elegans    | 165    | 87.4        | 38.64         | 14           | 1,00E-18 |
| g7861.t1   | TR | B4F769.1    | SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A-like protein 1 | Rattus norvegicus         | 910    | 251         | 37.21         | 99           | 2,00E-71 |
| g3679.t1   | TR | Q19269.2    | Zinc metalloproteinase nas-14                                                                       | Caenorhabditis elegans    | 503    | 6239        | 29.72         | 72           | 7,00E-22 |
| g15153.t1  | TR | P0DPW4.1    | U-scoloptoxin(01)-Er1a                                                                              | Ethmostigmus rubripes     | 119    | 68.9        | 39.76         | 46           | 2,00E-14 |
| g15484.t1  | TR | Q155U0.1    | Vacuolar protein sorting-associated protein 51 homolog                                              | Danio rerio               | 827    | 592         | 45.03         | 78           | 0        |
| g11960.t1  | TR | A0A4C1TJ24  | E3 ubiquitin protein ligase                                                                         | Eumeta variegata          | 917    | 1227        | 32.9          | 54.0         | 4.6e-150 |
| g8100.t1   | TR | Q2YDR3.1    | Inositol monophosphatase 3                                                                          | Danio rerio               | 341    | 173         | 36.04         | 22           | 3,00E-46 |
| g3428.t1   | TR | Q09510.1    | Myosin regulatory light chain                                                                       | Caenorhabditis elegans    | 172    | 199         | 56            | 91           | 9,00E-65 |


## Step 6. [Pfam](https://www.ebi.ac.uk/Tools/hmmer/) prediction. :pig2:
| Protein_id | Best blast hit | Relative organism                                    | Pham domains              | Probable localization (WoLF PSORT) | Probable localization (TargetP) |
|------------|----------------|------------------------------------------------------|---------------------------|------------------------------------|---------------------------------|
| g16368.t1  | 1e-05        | Penicillium italicum                                 | -                         | Nuclear                            | Other                           |
| g10514.t1  | 4e-13        | Pseudogymnoascus sp. VKM F-4246                      | -                         | Nuclear                            | Other                           |
| g14472.t1  | 0.0       | Denticeps clupeoides                                 | -                         | Nuclear                            | Other                           |
| g16318.t1  | 0.0            | Owenia fusiformis                                    | -                         | Nuclear                            | Other                           |
| g8312.t1   | 0.0            | Crassostrea gigas                                    | Clathrin                  | Nuclear                            | Other                           |
| g15484.t1  | 0.0        | Coptotermes formosanus                               | Vps51                     | Nuclear                            | Other                           |
| g7861.t1   | 2e-71        | Penicillium italicum                                 | SNF2-rel_dom, HARP        | Nuclear                            | Other                           |
| g11806.t1  | 3e-46        | Malassezia globosa (strain ATCC MYA-4612 / CBS 7966) | -                         | Nuclear                            | Other                           |
| g11960.t1  | 6e-98       | Eumeta variegata                                     | zf-C3HC4                  | Nuclear                            | Other                           |
| g8100.t1   | 3e-46        | Agrilus planipennis                                  | Inositol_P, MKLP1_Arf_bdg | Nuclear                            | Other                           |
| g5927.t1   | 1e-18        | Capitella teleta                                     | -                         | Nuclear                            | Other                           |
| g10513.t1  | 7e-22        | Malassezia globosa (strain ATCC MYA-4612 / CBS 7966) | -                         | Nuclear                            | Other                           |
