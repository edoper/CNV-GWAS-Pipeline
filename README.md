# CNV-GWAS-Pipeline

**CNV GWAS pipeline** describes a gold standard workflow to detect and analyze CNVs in large genome wide association studies (GWAS). We describe established CNV quality control procedures and case-control burden analysis for genotyping data. The present GitHub repository is associated to the "Genomic structural variants in nervous system disorders" chapter by Eduardo Pérez-Palma *et al* from the Neuromethods book series (2021).

## Dependencies
- [PLINK 1.9](https://www.cog-genomics.org/plink2)
- [KING](http://people.virginia.edu/~wc9c/KING/)
- [EIGENSTRAT](https://github.com/DReichLab/EIG)
- [PennCNV](http://penncnv.openbioinformatics.org/en/latest/)
- [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [R](https://www.r-project.org/)
- [Python](https://www.python.org/)
## Step 1. Quality control procedures.
For an extensive description of GWAS QC procedures please refer to [Anderson *et al*, 2010](https://www.nature.com/articles/nprot.2010.116) and [Clarke *et al*, 2011](https://www.nature.com/articles/nprot.2010.182). The format of genotyping data varies among SNP platforms. For the sake of simplicity, we asume that the input is on the standard PED and MAP file formats, which can be read by PLINK software.
We will transform the original data into PLINK BED files **(do not confuse this with UCSC BED format)** to facilitate the analyses of large datasets.
Remember to make backups of the original data. `plink --file <filename>` loads PED and MAP data. `plink --bfile <filename>` loads PLINK BED data.
```
plink --file mygwas --make-bed --out mygwas
```
With `(mygwas.map, mygwas.ped)` as the input, and `(mygwas.bed, mygwas.bim, mygwas.fam)` as the output.
### Genotype-level QCs
It is assumed that discordant sex information between the reported and genotyped sex in individuals is a sign of poor genotyping. We will check if there are discrepancies in sex and remove those individuals.
```
plink --bfile mygwas --check-sex --out mygwas
grep PROBLEM mygwas.sexcheck | awk '{print $1,$2}' > toremove.sexcheck.list
```
PLINK can filter data using multiple parameters. We will filter the data to remove samples with genotype call rate below 0.96 (`--mind 0.04`), SNVs with genotyping rate below 0.98 (`--geno 0.02`), minor allele frequency below 0.05 (`--maf 0.05`) and variants with a significant deviation from Hardy–Weinberg equilibrium (P < 0.001, `--hwe 0.001`). We will also remove the individuals that didn't pass the earlier test.
```
plink --bfile mygwas --remove toremove.sexcheck.list --mind 0.04 --geno 0.02 --maf 0.05 --hwe 0.001 --make-bed --out mygwas.genoQC
```
### Cohort-level QC
Cryptic relatedness and ancestry should be addressed to avoid spurious relationships in the analysis. KING will be used on the PLINK filtered output to identify relatedness in samples up to the second degree.
```
## UN = UNrelated
king -b mygwas.genoQC.bed --related –degree 2 --prefix mygwas.king
grep -v UN mygwas.king.kin0 | awk '{print $1,$2,$3}' > mygwas.related.list
```
If related individuals are found, check their genotype call rates (`--missing`) and remove the individual with less calling rates.

To assess ancestry, we merge the dataset with another known dataset with clearly-defined populations(`--merge` for non-binary format PLINK files, `--bmerge` for binary format PLINK files), and use PLINK to generate a Principal components analysis (PCA) of the different populations. [1000 genomes](https://www.internationalgenome.org/) is usually used to check for population stratification data.
```
## To download 1000 genomes data:
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
```
1000 genomes data is on VCF format, we need to convert it to PLINK binary format, filter it, and merge it with the original dataset.
Before merging, you need to make sure that the files are mergeable, for this there are three things that you need to keep in mind:
1) The reference genome needs to be the same in the dataset and the 1000 genomes data.
2) Resolve strand isues.
3) Remove the SNPs which still differ between the datasets.

The [PLINK merge](https://www.cog-genomics.org/plink/1.9/data#merge) manual has guidelines about what should be done if two PLINK binary datasets can not be merged.
```
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --biallelic-only strict --allow-no-sex --geno 0.02 --maf 0.05 --hwe 0.001 --make-bed --out 1000g.ALL
plink --bfile mygwas.genoQC --bmerge 1000g.ALL --make-bed --out mygwas.merged
plink --bfile mygwas.merged --pca --out mygwas.merged
```
The output (`mygwas.merged.eigenvec`) can be used to plot and identify outliers using R if the ancestry of every individual is appended as the last column of the PCA results.
```
library(readr)
mygwas_merged <- read_table2("mygwas.merged.eigenvec", col_names = FALSE)
library("ggplot2")
ggplot(data=testbed_merged_genoQC, aes(x=X3,y=X4,group=X23,color=X23)) + geom_point()
```
### Intensity-level QCs
TBD
## Step 2. CNV detection and downstream analysis.
### CNV detection
### Downstream analysis: Post-detection QCs.
Usually, CNV detection is restricted to the autosomal chromosomes as the intensity analysis of probes mapping within the X and Y chromosome is less reliable.
We can remove those chromosomes from the analysis with `grep -v "chrX"` and `grep -v "chrY"`.
We also need to filter CNVs with less than 20 probes, a genomic area below 20 kb, and a marker density (probes/genomic area) below 0.0001. We can filter the data using Python, as shown in the following example:
```
## Save as cnvfilter.py
## Usage: python cnvfilter.py <cnvlist>
from sys import argv

with open(argv[1]) as f1:
	for i in f1:
		probes = float(i.strip().split()[1].split("=")[1])
		distance = float(i.strip().split()[2].split("=")[1].replace(',', ''))
		if probes > 20:
			if distance > 20000:
				if (probes/distance) > 0.0001:
					print(i.strip())
```
After filtering, we strongly recommend to implement a quality score calculation for each CNV following the methods proposed by [Macé et al. 2016](https://academic.oup.com/bioinformatics/article/32/21/3298/2415363), in which various CNV metrics are combined to estimate the probability of a called CNV to be a consensus call. 
A similar script can be used to transform the data into UCSC BED format.
```
awk '{print $1}' | awk -F":|-" '{print $1,$2,$3}' OFS="\t"
```
False positive CNVs calls tend to fall within highly repetitive regions such as the centromere, telomer and immunoglobulins regions and calls overlapping the boundaries of repetitive regions should be removed. This step can be done with bedtools, if the unwanted regions are in UCSC BED format.
```
intersectBed -a all-cnv.bed -b filterregions.bed -v > all-cnv.filtered.bed
```
### Downstream analysis: Annotation
The annotation of the QC-passing CNVs is essential to extract significant biological knowledge from a case-control cohort. The filtered `all-cnv.filtered.bed` UCSC BED file can be directly uploaded to Ensembl's variant effect predictor [VEP](https://www.ensembl.org/Tools/VEP). However, for larger annotation procedures local VEP installation is recommended. [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/) or PennCNV can also be used to annotate the dataset.
## Step 3. Burden analysis.
### Concept
### Sample-wise input.
### Binary CNV predictor
### Logistic regression
