# CNV-GWAS-Pipeline

**CNV GWAS pipeline** describes a gold standard workflow to detect and analyze CNVs from large genome wide association studies (GWAS). We describe established CNV quality control procedures and case-control burden analysis for genotyping data. The present GitHub repository is associated to the ["Genomic structural variants in nervous system disorders"](https://not.ready.yet) chapter by Eduardo Pérez-Palma *et al* from the Neuromethods book series (2021). *Disclaimer: The data presented here is for educational purpuses only and does not represent real results or association*.   

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
For a detailed description of GWAS quality control (QC) procedures please refer to [Anderson *et al*, 2010](https://www.nature.com/articles/nprot.2010.116) and [Clarke *et al*, 2011](https://www.nature.com/articles/nprot.2010.182). The format of genotyping data varies among SNP platforms. For simplicity, we asume that the input is on the standard PED and MAP file formats, which can be read by PLINK software.
We will transform the original data into PLINK BED files **(do not confuse this with UCSC BED format)** to facilitate the analyses of large datasets.
Remember to make backups of the original data. The command `plink --file <filename>` loads PED and MAP data. The command `plink --bfile <filename>` loads PLINK BED data.
```
plink --file mygwas --make-bed --out mygwas
```
With `(mygwas.map, mygwas.ped)` as the input, and `(mygwas.bed, mygwas.bim, mygwas.fam)` as the output.
### Genotype-level QCs
It is assumed that discordant sex information between the reported and genotyped sex in individuals is a sign of poor genotyping. The following command will check if there are discrepancies in sex and remove those individuals.
```
plink --bfile mygwas --check-sex --out mygwas
grep PROBLEM mygwas.sexcheck | awk '{print $1,$2}' > toremove.sexcheck.list
```
PLINK can filter data using multiple parameters. We will filter the data to remove samples with genotype call rate below 0.96 (`--mind 0.04`), SNVs with genotyping rate below 0.98 (`--geno 0.02`), minor allele frequency below 0.05 (`--maf 0.05`) and variants with a significant deviation from Hardy–Weinberg equilibrium (P < 0.001, `--hwe 0.001`). We will also remove the individuals that didn't pass the earlier test. Such a command states:
```
plink --bfile mygwas --remove toremove.sexcheck.list --mind 0.04 --geno 0.02 --maf 0.05 --hwe 0.001 --make-bed --out mygwas.genoQC
```
### Cohort-level QC
Cryptic relatedness and ancestry should be addressed to avoid spurious relationships in the analysis. KING will be used on the PLINK filtered output to identify relatedness in samples up to the second degree. Below commands will identify relted individuals in mygwas data.  
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
You need to make sure that the files are mergeable, for this there are three things that you need to keep in mind:
1) The reference genome needs to be the same in your dataset and in the 1000 genomes data.
2) Resolve strand isues.
3) Remove the SNPs which still differ between the datasets.

The [PLINK merge](https://www.cog-genomics.org/plink/1.9/data#merge) manual has guidelines about what should be done if two PLINK binary datasets can not be merged. The below command will merge mygwas with 1000 genomes data and carry out a PCA analysis. 
```
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --biallelic-only strict --allow-no-sex --geno 0.02 --maf 0.05 --hwe 0.001 --make-bed --out 1000g.ALL
plink --bfile mygwas.genoQC --bmerge 1000g.ALL --make-bed --out mygwas.merged
plink --bfile mygwas.merged --pca --out mygwas.merged
```
The output (`mygwas.merged.eigenvec`) can be used to plot and identify outliers using R if the ancestry of every individual is appended as the last column of the PCA results. The command below in R will plot mygwas individuals alongside 1000 genomes samples. 
```
library(readr)
mygwas_merged <- read_table2("mygwas.merged.eigenvec", col_names = FALSE)
library("ggplot2")
ggplot(data=testbed_merged_genoQC, aes(x=X3,y=X4,group=X23,color=X23)) + geom_point()
```
After filtering outliers, we will require a table with the individual ID, phenotype, sex, and the first three components from the PCA.
The following Python script can be used to generate this table, with the PLINK fam file and the PLINK eigenvec file as input.
```
## Save as pcatablebuilder.py
## Usage: python pcatablebuilder.py <fam indlist> <PCA results.eigenvec>

from sys import argv
print("SampleID	Response	Sex	PC1	PC2	PC3")

with open(argv[1]) as f1:
	for i in f1:
		splitfirst = i.strip().split()
		
		with open(argv[2]) as f2:
			for j in f2:
				splitsecond = j.strip().split()
				if splitfirst[1] == splitsecond[1]:
					joined = "	".join(splitfirst[1],splitfirst[-1],splitfirst[-2],splitsecond[2],splitsecond[3],splitsecond[4])
					print(joined)

```
```
python pcatablebuilder.py mygwas.fam mygwas.merged.eigenvec > mysamples
```

## Step 2. CNV detection and downstream analysis.
### CNV detection
For CNV detection, we strongly recommend to follow the pipeline proposed by [Macé et al. 2016](https://academic.oup.com/bioinformatics/article/32/21/3298/2415363), which uses PennCNV to detect CNVs and R scripts for filtering. The pipeline describes how to install PennCNV, R, and the required R libraries to properly execute the R scripts. The pipeline requires a configuration file.
```
## config_example.txt
pennCNVpath:    /path/to/pennCNV/
HMMpath:	/path/to/pennCNV/lib/hhall.hmm
HMMcreate:	0
PFB:	/path/to/input/pfb/
CompilePFB:	1
GCmod:	my/GC_model/path
UseGCmod:	0
InputData:	1
DATA:	/path/to/input/
OUTPUT:	/path/to/output/results
FormatedPath:	/path/to/input/formated
Chromosome:	1-22
CNVcall:	1
Cleancall:	1
format:	1
CreateRfile:	1
AssoData:	1
NbCores:	16
PhenoPath:	/path/to/input/phenotype/phenotype.txt
Phenotype:	phenoName
```
```
./CNV_detection.sh config_example.txt
```
The PennCNV Pipeline User Guide included in the pipeline describes the complete output, here we highlight two files:

- mygwas.clean.rawcnv contains the ALL detected CNVs
- mygwas.good.cnv contains the filtered subset of CNVs in PennCNV format and is the main input for downstream analysis.

Here an example of the PennCNV calls output:
Genomic coordinates | N SNPs | Lenght | Copy number (Del=0/1;Dup=3/4) | Sample | Begin | End
--- | --- | --- | --- | --- | --- | ---
chr1:25593128-25611452 | numsnp=14 | length=18,325 | state1,cn=0 | NA19222 | startsnp=CN_482242 | endsnp=CN_020771
chr1:72771143-72811148 | numsnp=44 | length=40,006 | state1,cn=0 | NA21596 | startsnp=CN_517829 | endsnp=CN_519942
chr1:152555795-152586594 | numsnp=33 | length=30,800 | state1,cn=0 | NA20787 | startsnp=CN_452251 | endsnp=CN_453516

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
```
python cnvfilter.py mygwas.good.cnv >mygwas.good.filtered.cnv
```
Next to remove custom regions we need a bed version of our CNV calls. An awk one-liner can be used to transform the data into UCSC BED format.
```
less mygwas.good.filtered.cnv | awk '{print $1"\t"$5}' | perl -pe 's/:|-/\t/g' >mygwas.good.filtered.cnv.bed
perl -pi -e 's/\Achr//g' mygwas.good.filtered.cnv.bed
sort -V -k1,2 mygwas.good.filtered.cnv.bed >mygwas.good.filtered.cnv.sorted.bed
```
False positive CNVs calls tend to fall within highly repetitive regions such as the centromere, telomer and immunoglobulins regions. Calls overlapping the boundaries of repetitive regions should be removed. Here, we will exclude [repetitive regions](https://github.com/dellytools/delly/blob/master/excludeTemplates/human.hg19.excl.tsv) with bedtools (see repetitive-regions bed file). A single line will exclude undesired CNVs:
```
intersectBed -a mygwas.good.filtered.cnv.sorted.bed -b repetitive-regions.bed -v >mycnvs.final.bed
```
This way we obtain our final set of CNV calls mycnv-final.bed. Note that this file contains one cnv per row. If same cnv was called on three samples, three rows will show the same cnv.  See example below:

mycnv-final.bed first three rows: 

Chr | Start | End | Sample
--- | --- | --- | ---
1 | 61735 | 235938 | NA07435
1 | 61735 | 235938 | NA11839
1 | 61735 | 235938 | NA12348

### Downstream analysis: Annotation
The annotation of the QC-passing CNVs is essential to extract significant biological knowledge from a case-control cohort. For simplicity in this example we will annotate all canonical RefSeq genes (see refseq.genes.bed file).
```
bedtools intersect -a mycnvs.final.bed -b refseq.genes.bed -wa -wb | awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8}' >mycnvs.intersected
cat mycnvs.intersected | bedtools groupby -g 1,2,3,4 -c 5 -o count,collapse >mycnvs.annotated
```
The expected output should look like the table below: 

Chr | Start | End | Sample | N genes | Gene Names
--- | --- | --- | --- | --- | --- | 
1 | 14879491 | 14951408 | NA20760 | 1 | KAZN
1 | 15718470 | 15789733 | NA21304 | 3 | CTRC,EFHD2,CELA2A
1 | 15894607 | 16000741 | NA18534 | 4 | RSC1A1,AGMAT,DNAJC16,DDI2

For further annotations, you can directly upload the filtered `all-cnv.filtered.bed` BED file to the Ensembl's variant effect predictor [VEP](https://www.ensembl.org/Tools/VEP) to annotate all relevant biological features. However, for larger annotation procedures local VEP installation is recommended. [ANNOVAR](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/) or PennCNV can also be used to annotate the dataset.

## Step 3. Burden analysis.

### Regions of Interest (or genes of interest)
CNV burden analysis is a hypothesis-driven approach that requires the definition of at least one region of interest. A region of interest can be any genomic interval defined by the user, and usually takes the form of a gene set. Gene sets should be meaningful to the disease to extract valuable conclusions. There is no restriction in the type or number of regions that can be tested. 

To annotate regions of interest you might take advantage of the refseq.genes.bed file to generate a regions.of.interest.bed file. For this example we will as "regions of interest" a list of 279 coding genes (gene.set.list) known to be associated to developmental disorders. Next you can interrogate is any CNVs are overlapping these regions and generate a mysamples.w.predictor list.  
```
awk -F"\t" 'BEGIN{OFS="\t"} FNR==NR{p[$4]=$0;next}{print p[$1]}' refseq.genes.bed gene.set.list | sort -V -k1,2 >regions.of.interest.bed
bedtools intersect -a mycnvs.final.bed -b regions.of.interest.bed -wa | awk '{print $4"\t1"}' | sort | uniq >mysamples.w.predictor
```

### Input.
To carry out a burden analysis you need three variables for every sample included in the gwas:
1. RESPONSE: binary phenotype (cases=1; controls=0) in all PostQC gwas samples regardless of their cnv state.
2. PREDICTOR: binary precdictor, in this case if a sample has or not a CNV overlapping a region of interest. we will use mysamples.w.predictor to create the PREDICTOR field.  
3. COVARIABLES (optional): Here you can inlcude sex or principal components.

RESPONSE and COVARIABLES were extracted from STEP1 (mysamples file). Here, the PREDICTOR needs to be construted. we simply add a mock negative PREDICTOR column and then map the true value with the file mysamples.w.predictor generated in the previous step.
```
awk -F"\t" 'BEGIN{OFS="\t"} FNR==NR{p[$1]=$2;next}{print $0, p[$1]}' mysamples.w.predictor mysamples >sample.wise.input
perl -pi -e 's/\t\n/\t0\n/g' sample.wise.input
```
The final sample wise input shoudl look like 

Sample | RESPONSE | Sex | PC1 | PC2 | PC3 | PREDICTOR
--- | --- | --- | --- | --- | --- | ---
NA06984 | 0 | 1 | 158123 | 391066 | 139324 | 1
NA06989 | 1 | 1 | 158886 | 387969 | 146544 | 1
NA12335 | 0 | 1 | 159503 | 386152 | 140885 | 0

### Logistic regression
Having sample.wise.input with PREDICTOR annotatated the user can proceed with the logistic regression of the burden analysis. The following commands will test for association between the RESPONSE (i.e. phenotype) and the PREDICTOR (i.e. having a CNV in a region on interest). 

```
#load libraries
library(aod)
library(Rcpp)
library(MASS)

#load data
data.cnv <- read.table("sample.wise.input", header=T)

#model 
model.burden <- glm(RESPONSE ~ PREDICTOR + Sex + PC1 +PC2 +PC3, data=data.cnv, family=binomial)

#data extraction
sum.mod <- summary(model.burden)
sum.tab <-sum.mod$coefficients
p.mod<-with(model.burden, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
or.ci.mod<-exp(cbind(OR = coef(model.burden), confint(model.burden)))

#results
main.results<-cbind(p.mod, or.ci.mod)
main.results
```
Expected results should look like the tabe below: 

Variables | p.model | OR | 2.5 % | 97.5 %
--- | --- | --- | --- | --- 
(Intercept) | 5.339691e-19 | 0.6600314 | 0.5455037 | 0.7971419
PREDICTOR | 5.339691e-19 | 2.9084473 | 2.3352957 | 3.6305632
Sex | 5.339691e-19 | 0.9552320 | 0.7683386 | 1.1873764
PC1 | 5.339691e-19 | 1.0000002 | 0.9999998 | 1.0000007
PC2 | 5.339691e-19 | 1.0000001 | 0.9999998 | 1.0000004
PC3 | 5.339691e-19 | 0.9999999 | 0.9999996 | 1.0000002

If you are able to obtain the above table, Congratulations! you have just carried out a CNV burden analysis based on GWAS data. In our example, there is significant association between the phenotype of our cases and CNVs overlapping genes associated to developmental disorders. In other words, the burden of patients with CNVs overlapping genes associated to developmental disorders is significantly higher than the one observed in controls.

Thank you!
