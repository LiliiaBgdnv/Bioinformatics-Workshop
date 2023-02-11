# Progect 5.Human genome analysis.

Create new environment:

```ruby
python3.11 -m venv pl
source progect5/bin/pl
```

Downloaded a dataset with data and [Plink](https://www.cog-genomics.org/plink/)

```ruby
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip
sudo chmod +x plink
```

Remove all SNPs corresponding to deletions and insertions to make the file compatible with annotation tools.

```ruby
./plink --23file SNP_raw_v4_Full_20170514175358.txt --recode vcf --out snps_clean --output-chr MT --snps-only 'just-acgt'
```

## Identification of haplogroups [Morleydna](https://ytree.morleydna.com/extractFromAutosomal)
[Result](file:///C:/Users/iskys/AppData/Local/Temp/Temp1_processedData.zip/processedData/MorleyDNA.com_output/MorleyDNA.com%20Y-SNP%20Subclade%20Predictor.html)

## Annotation - sex and eye colour
In our data there is Y chromosome, so it can be assumed that the genome belongs to a man.

We determined eye color using data from the [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694299/)

![image](https://user-images.githubusercontent.com/109213422/218263729-a117a5e3-dd50-49fb-9d55-166dd04f50b9.png)

rs12913832 **AG** => not blue; rs16891982 **CG**; rs6119471 not in test; rs12203592 **CT**; rs16891982 **CG**

We couldn't tell exactly, so we assumed the eyes were **green-brown**.

## Annotation of all SNPs, selection of clinically relevant ones using VEP ([Variant Effect Predictor](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP)).

Fildering result:
```ruby
awk '($32!="-") ' vep_result.txt | grep risk_factor | cut -f 1-3 | sort | uniq > vep_filter.txt
```


| Uploaded variant | Location | Allele  | Meaning                                                                                                                          |
|------------|------------------------|------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| rs1024611  | 17:32579788-32579788   | G          | discussed in the 23andMe blog as being relevant to HIV                                                                                                                         |
| rs1049296  | 3:133494354-133494354  | T          | Involved in the transport of iron, transferrin and its subtypes have been linked at times to various conditions, perhaps most notably Alzheimer's disease.                     |
| rs10757274 | 9:22096055-22096055    | G          | ~1.3x increased risk for heart disease                                                                                                                                         |
| rs1169288  | 12:121416650-121416650 | C          |  The genetic susceptibility to type 2 diabetes ???                                                                                                                             |
| rs12150220 | 17:5485367-5485367     | T          | slightly increased risk for several autoimmune diseases                                                                                                                        |
| rs13266634 | 8:118184783-118184783  | T          | increased risk for type-2 diabetes                                                                                                                                             |
| rs1801275  | 16:27374400-27374400   | G          | 1.4x higher risk for meningiomas                                                                                                                                               |
| rs2004640  | 7:128578301-128578301  | T          | 1.4x increased risk for SLE (волчанка)                                                                                                                                         |
| rs2073658  | 1:161010762-161010762  | T          | is associated with a modestly increased risk to develop type 2 diabetes                                                                                                        |
| rs2241880  | 2:234183368-234183368  | G          | strongly associated with ileal Crohn's disease                                                                                                                                 |
| rs231775   | 2:204732714-204732714  | G          | 2.3x risk of Hashimoto's thyroiditis, 1.47x risk of Graves' disease                                                                                                            |
| rs4402960  | 3:185511687-185511687  | T          | 1.2x increased risk for type-2 diabetes, ~1x risk for gestational diabetes                                                                                                     |
| rs4961     | 4:2906707-2906707      | T          | 1.8x increased risk for high blood pressure                                                                                                                                    |
| rs4977574  | 9:22098574-22098574    | G          | Most studies find a somewhat elevated (~1.5x) risk for myocardial infarction                                                                                                   |
| rs5186     | 3:148459988-148459988  | C          | ~1.4x increased risk of hypertension                                                                                                                                           |
| rs61747071 | 16:53720436-53720436   | T          | uncertain-significance,                                                                                                                                                        |
| rs6265     | 11:27679916-27679916   | T          | Slightly increased risk for ADHD or depression; somewhat quicker mental decline in Alzheimer patients                                                                          |
| rs6280     | 3:113890815-113890815  | T - normal | rs6280, also known as Ser9Gly, is a SNP in the dopamine receptor D3 DRD3 gene, those who were rs6280(C;C) homozygotes had greater positive symptom remission for schizophrenia |
| rs699      | 1:230845794-230845794  | G          | increased risk of hypertension                                                                                                                                                 |
| rs763110   | 1:172627498-172627498  | T          | ~0.80x reduced cancer risk                                                                                                                                                     |
| rs7794745  | 7:146489606-146489606  | T          | slightly increased risk for autism                                                                                                                                             |
| rs909253   | 6:31540313-31540313    | G          | Myocardial infarction Psoriatic arthritis                                                                                                                                      |
