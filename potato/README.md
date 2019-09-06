---
title: "Building a genetic map of an autotetraploid population using MAPpoly"
author: "Marcelo Mollinari, Gabriel Gesteira, Guilhereme Pereira, A Augusto F Garcia and Zhao-Bang Zeng"
date: '2019-09-06'
output:
  html_document:
    highlight: pygments
    keep_md: yes
    theme: flatly
    toc: yes
    toc_depth: '3'
    toc_float:
      collapsed: no
  md_document:
    variant: markdown_github
  pdf_document:
    highlight: pygments
    toc: yes
    toc_depth: '2'
  word_document:
    toc: yes
    toc_depth: '2'
linestretch: 1.2
bibliography: biblio.bib
---



# Introduction

`mappoly` (v. 0.1.0) is an R package to construct genetic maps in autopolyploids with even ploidy levels. In its current version, `mappoly` can handle ploidy levels up to 8 when using hidden Markov models (HMM), and up to 12 when using the two-point simplification. All the two-point based functions are fast enough to run on standard computers. However, we strongly recommend using high-performance computation for HMM-based analysis, especially for ploidy levels higher than 4. 

Here we assume that the genotypic data is available and in the format required by `mappoly`. The primary purpose of this tutorial is to show some functions available in `mappoly` and how to use them sequentially to construct a genetic map. TThe derivation of the HMM used in `mappoly` can be found in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378).

`mappoly` is not yet available in CRAN, but you can install it from Git Hub. Within R, you need to install the package `devtools`:

```R
install.packages("devtools")
```
To install `mappoly` from Git Hub use

```R
devtools::install_github("mmollina/mappoly")
```

# Loading `mappoly`

To load `mappoly`, simply type 


```r
library(mappoly)
```

```
## Registered S3 methods overwritten by 'car':
##   method                          from
##   influence.merMod                lme4
##   cooks.distance.influence.merMod lme4
##   dfbeta.influence.merMod         lme4
##   dfbetas.influence.merMod        lme4
```

## Loading datasets

In its current version, MAPpoly can handle three different types of datasets:

- CSV files
- MAPpoly files
  - Dosage based
  - Probability based
- VCF files (beta)

Please notice that both CSV and MAPpoly datasets are sensible to common errors, such as additional spaces, commas and wrong encoding (non-UTF-8). If you have trouble, please double check your files before submitting an issue. Detailed steps of all supported files are given on the sections below. Support for Illumina files is coming soon.

### Reading CSV files

The preparation of a CSV file for MAPpoly is quite straightforward. It can be done in Microsoft Excel or any other spreadsheet software of your preference.  In this file, each line comprehends a marker and each column comprehends information about the marker. In its current version, MAPpoly can handle .csv files with allelic dosage data. 

The first line of the CSV file should contain headers for all columns. The first five columns should contain the following information: marker name, the dosage of both parents, a sequence number (e.g., a chromosome number, if available) and a sequence position (e.g., the marker position within the chromosome, if available). In addition to these five headers, you should include the name of all individual in the population. From the second line onwards, all columns should contain its values, including allelic dosages for all individuals. Missing or absent values should be represented by NA.

NOTE: If genomic information is not available, the 'sequence' and 'sequence position' columns should be filled with NA's.

Example:

<div class="figure" style="text-align: center">
<img src="image/csv.png" alt="Figure 1: Example of CSV data set" width="90%" />
<p class="caption">Figure 1: Example of CSV data set</p>
</div>

**Important note: avoid spaces in .csv files.** As mentioned above, please double check your datasets for extra spaces, commas, dots and encoding. Your .csv file should be encoded using UTF-8.

You can read csv files with the `read_geno_csv` function:


```r
solcap.file <- system.file('extdata', 'tetra_solcap.csv', package = 'mappoly')
dat.dose.csv <- read_geno_csv(file.in  = solcap.file, ploidy = 4)
```

```
## Reading the following data:
##     Ploidy level: 4
##     No. individuals:  160
##     No. markers:  4017
##     No. informative markers:  4017 (100%)
##     This dataset contains sequence information.
##     ...
##     Done with reading.
##     Filtering non-conforming markers.
##     ...
##     Done with filtering.
```


In addition to the .csv file path, you should indicate the ploidy level using the `ploidy` argument. This function automatically excludes monomorphic markers, keeping only informative ones. It also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default: TRUE), which excludes non-expected genotypes under these assumptions.

### Reading MAPpoly files

Besides CSV and VCF files, MAPpoly can also handle two more dataset types that follow the same format: a genotype-based file (wth allele dosages) (1) and a probability-based file (2). Both are pure text files with the same header, but different genotype table formats.

For both the header should contain ploidy level, number of individuals (nind), number of markers (nmrk), marker names (mrknames), individual names (indnames), allele dosages for parent 1 (dosageP), allele dosages for parent 2 (dosageQ), sequence/chromosome information (seq), position of each marker (seqpos), number of phenotipic traits (nphen) and the phenotypic data (pheno) if available. The header should be organized according to this example:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
```

For more information about MAPpoly file format, please see `?read_geno` and `?read_geno_dist` documentation from MAPpoly package.

#### Using `read_geno`

The header should be followed by a table containing the genotypes (allele dosages) for each marker (rows) and for each individual (columns), as follows:

|          | Individual 1 | Individual 2 | Individual 3 |
|----------|:------------:|:------------:|:------------:|
| Marker 1 | 1            | 0            | 0            |
| Marker 2 | 3            | 0            | 2            |
| Marker 3 | 1            | 0            | 0            |
| Marker 4 | 1            | 0            | 0            |
| Marker 5 | 3            | 4            | 4            |

The final file should look like the example below:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
1 0 0
3 0 2
1 0 0
1 0 0
3 4 4
```

Then, use the `read_geno` function to read your file:


```r
solcap.geno.file <- system.file('extdata', 'tetra_solcap_geno', package = 'mappoly')
dat.dose.mpl <- read_geno(file.in  = solcap.geno.file)
```

```
## Reading the following data:
##     Ploidy level: 4
##     No. individuals:  160
##     No. markers:  4017
##     No. informative markers:  4017 (100%)
##     This dataset contains sequence information.
##     ...
## 
##     Done with reading.
##     Filtering non-conforming markers.
##     ...
##     Performing chi-square test.
##     ...
##     Done.
```

```r
dat.dose.mpl
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             4017 
##     Missing data:                            2.93%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  434
##      1  1  309
##      1  2  256
##      1  3  125
##      1  4   12
##      2  0  171
##      2  1  212
##      2  2  337
##      2  3  196
##      2  4   50
##      3  0   41
##      3  1  110
##      3  2  202
##      3  3  203
##      3  4  334
##      4  1    9
##      4  2   86
##      4  3  309
```

This function automatically excludes monomorphic markers, keeping only informative ones. It also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default: TRUE), which excludes non-expected genotypes under these assumptions. You can also define the p-value threshold used in the segregation test with the `thresh.line` parameter. 

#### Using `read_geno_dist`

Following the same header (described before) should be a table containing the probability distribution for each combination of marker $\times$ individual. Each line of this table represents the combination of one marker with one individual, and its respective probabilities of having each possible allele dosage. The first two columns represent the marker and the individual, respectively, and the remaining elements represent the probability associated with each one of the possible dosages, as follows:

| Marker | Individual | $p(d=0)$ | $p(d=1)$ | $p(d=2)$ | $p(d=3)$ | $p(d=4)$ |
|--------|:----------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| M1     | Ind1       | 0.5      | 0.5      | 0.0      | 0.0      | 0.0      |
| M2     | Ind1       | 0.0      | 1.0      | 0.0      | 0.0      | 0.0      |
| M3     | Ind1       | 0.3      | 0.7      | 0.0      | 0.0      | 0.0      |
| M4     | Ind1       | 0.5      | 0.5      | 0.0      | 0.0      | 0.0      |
| M5     | Ind1       | 0.0      | 0.0      | 0.0      | 0.9      | 0.1      |
| M1     | Ind2       | 1.0      | 0.0      | 0.0      | 0.0      | 0.0      |
| M2     | Ind2       | 0.2      | 0.5      | 0.3      | 0.0      | 0.0      |
| M3     | Ind2       | 0.9      | 0.1      | 0.0      | 0.0      | 0.0      |
| M4     | Ind2       | 0.9      | 0.1      | 0.0      | 0.0      | 0.0      |
| M5     | Ind2       | 0.0      | 0.0      | 0.0      | 0.2      | 0.8      |
| M1     | Ind3       | 0.2      | 0.8      | 0.0      | 0.0      | 0.0      |
| M2     | Ind3       | 0.4      | 0.6      | 0.0      | 0.0      | 0.0      |
| M3     | Ind3       | 1.0      | 0.0      | 0.0      | 0.0      | 0.0      |
| M4     | Ind3       | 0.0      | 0.1      | 0.9      | 0.0      | 0.0      |
| M5     | Ind3       | 0.1      | 0.9      | 0.0      | 0.0      | 0.0      |

Please notice that each marker $\times$ individual combination have $m+1$ associated probabilities, being $m$ the ploidy level and $m+1$ the number of possible allele dosages. The final file (header + table) should look like this:

```
ploidy 4
nind 3
nmrk 5
mrknames M1 M2 M3 M4 M5
indnames Ind1 Ind2 Ind3
dosageP 0 2 0 0 3
dosageQ 1 2 1 1 3
seq 1 1 2 2 3
seqpos 100 200 50 150 80
nphen 0
pheno-----------------------
geno------------------------
M1 Ind1 0.5 0.5 0.0 0.0 0.0
M2 Ind1 0.0 1.0 0.0 0.0 0.0
M3 Ind1 0.3 0.7 0.0 0.0 0.0
M4 Ind1 0.5 0.5 0.0 0.0 0.0
M5 Ind1 0.0 0.0 0.0 0.9 0.1
M1 Ind2 1.0 0.0 0.0 0.0 0.0
M2 Ind2 0.2 0.5 0.3 0.0 0.0
M3 Ind2 0.9 0.1 0.0 0.0 0.0
M4 Ind2 0.9 0.1 0.0 0.0 0.0
M5 Ind2 0.0 0.0 0.0 0.2 0.8
M1 Ind3 0.2 0.8 0.0 0.0 0.0
M2 Ind3 0.4 0.6 0.0 0.0 0.0
M3 Ind3 1.0 0.0 0.0 0.0 0.0
M4 Ind3 0.0 0.1 0.9 0.0 0.0
M5 Ind3 0.1 0.9 0.0 0.0 0.0
```

Once the file is ready, use the function `read_geno_dist` to read it:


```r
solcap.file <- system.file('extdata', 'tetra_solcap_geno_dist.bz2', package = 'mappoly')
dat.dist.mpl <- read_geno_dist(file.in  = solcap.file, prob.thres = 0.95)
```

```
## Reading the following data:
##     Ploidy level: 4
##     No. individuals:  160
##     No. markers:  4017
##     No. informative markers:  4017 (100%)
##     This dataset contains sequence information.
##     ...
##     Done with reading.
##     Filtering non-conforming markers.
##     ...
##     Performing chi-square test.
##     ...
##     Done.
```

```r
dat.dist.mpl
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             4017 
##     Missing data under 0.95 prob. threshold: 2.93%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  434
##      1  1  309
##      1  2  256
##      1  3  125
##      1  4   12
##      2  0  171
##      2  1  212
##      2  2  337
##      2  3  196
##      2  4   50
##      3  0   41
##      3  1  110
##      3  2  202
##      3  3  203
##      3  4  334
##      4  1    9
##      4  2   86
##      4  3  309
```


**Important note:** as this type of file contains the probability distribution of the genotypes, it may take longer to read. 

This function automatically excludes monomorphic markers, keeping only informative ones. You can define the minimum probability necessary to call a dosage using the `prob.thres` argument. If the higher probability for a maker $\times$ individual passes this threshold, then its associated dosage is used. However, if none of the probabilities reach this threshold, then its dosage is considered missing (NA).
This function also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default: TRUE), which excludes non-expected genotypes under these assumptions. <!-- You can also define the p-value threshold used in the segregation test with the `thresh.line` parameter. -->

### Reading VCF files

VCF files are less sensible to errors, because they are usually produced by automated SNP calling pipelines and less susceptible to user edition. MAPpoly can also handle VCF files (version 4.0 and higher) produced by the most common softwares such as TASSEL, GATK, Stacks and many others. However, few available softwares can handle poliploidy and estimate genotypes (allele dosages) correctly. Due to this, many other softwares are dedicated to correctly estimate the allele dosages considering ploidy. Briefly, these softwares model the ratio between allele read counts for each marker $\times$ individual combination, and determines which is the most probable allele dosage given the observed ratio and other *a priori* information. Examples of these softwares are [SuperMASSA](http://statgen.esalq.usp.br/SuperMASSA/), [fitTetra](https://cran.r-project.org/web/packages/fitTetra/index.html), [ClusterCall](https://potatobreeding.cals.wisc.edu/wp-content/uploads/sites/161/2017/08/ClusterCall_Download.zip), [updog](https://github.com/dcgerard/updog), [SNPready](https://github.com/cran/snpReady) and many others. After allele dosage estimation, your VCF file should contain values like **1/1/1/0** (for an autotetraploid) rather than **1/0**. Since MAPpoly still needs allele dosages (or their probabilities) to build genetic maps, we strongly recommend that you use one of these softwares to estimate allele dosages before building the map. We are working on an automated integration of these softwares with MAPpoly.

Once you have your VCF file, read it with `read_vcf` function:


```r
tetra_vcf_file = system.file('extdata', 'tetra_example.vcf.gz', package = 'mappoly')
dat.dose.vcf = read_vcf(file.in = tetra_vcf_file, parent.1 = 'PARENT1', parent.2 = 'PARENT2')
```

```
## Reading data...
## Scanning file to determine attributes.
## File attributes:
##   meta lines: 8
##   header_line: 9
##   variant count: 9345
##   column count: 248
## Meta line 8 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
##   Character matrix gt rows: 9345
##   Character matrix gt cols: 248
##   skip: 0
##   nrows: 9345
##   row_num: 0
## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant: 9345
## All variants processed
## Processing genotypes...Done!
## Done!
## Read the following data:
##     Ploidy level: 4
##     No. individuals:  237
##     No. markers:  9345
##     No. informative markers:  4610 (49.3%)
##     This dataset contains sequence information.
##     ...
##     Done with reading.
##     Filtering non-conforming markers.
##     ...
##     Performing chi-square test.
##     ...
##     Done.
```

```r
dat.dose.vcf
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         237 
##     No. markers:                             4610 
##     Missing data:                            15.86%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1 1132
##      0  2  583
##      0  3   48
##      1  0  338
##      1  1 1148
##      1  2  582
##      1  3   51
##      1  4    6
##      2  0  171
##      2  1  233
##      2  2  198
##      2  3   12
##      2  4    1
##      3  0   80
##      3  1   19
##      3  2    1
##      3  4    1
##      4  1    5
##      4  2    1
```


<!-- **Important note:** function `read_vcf` is under development, and still takes a long time to read -->
Besides the path to your VCF file, you should indicate `parent.1` and `parent.2` names. Please notice that their names should be exactly the same strings that appear in your VCF file. The ploidy level will be automatically detected, but you may indicate it using the `ploidy` argument (optional) to let the function check for possible errors. This function also performs chi-square tests for all markers, considering the expected segregation patterns under Mendelian inheritance, random chromosome pairing and no double reduction. You can optionally use the `filter.non.conforming` logical argument (default: TRUE), which excludes non-expected genotypes under these assumptions. The p-value threshold used the segregation test can be defined by the `thresh.line` argument. <!-- Arguments `update.prob` and `output` are coming soon. -->

The example .vcf dataset from _Urochloa decumbens_ [Ferreira2019](https://www.frontiersin.org/articles/10.3389/fpls.2019.00092/full) was obtained from raw sequencing data available at [NCBI](https://www.ncbi.nlm.nih.gov/sra/SRP148665), using both [TASSEL4Poly](https://github.com/guilherme-pereira/tassel4-poly) and [VCF2SM pipeline](https://github.com/guilherme-pereira/vcf2sm).

### Combining multiple datasets

Function `merge_data` is coming soon.

## Exploratory Analysis 

For didatic purposes, we will keep using the tetraploid potato array data (loaded using the examples above). We will construct a genetic map of the B2721 population, which is a cross between two tetraploid potato varieties: Atlantic and B1829-5. The population comprises 160 offsprings genotyped with the SolCAP Infinium 8303 potato array. The dataset also contains the genomic order of the SNPs from the _Solanum tuberosum_ genome version 4.03. The genotype calling was performed using fitTetra R package using [this pipeline](https://github.com/mmollina/Autopolyploid_Linkage/blob/master/src/solcap_map_construction/snp_calling/genotype_calling_public_data_fittetra.R). Another option would be to use ClusterCall and [this pipeline](https://mmollina.github.io/tutorials/solcap/solcap_example.html).

Once the data is loaded, you can explore the dataset using the `print` function:


```r
print(dat.dose.mpl, detailed = TRUE)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             4017 
##     Missing data:                            2.93%
## 
##     ----------
##     No. markers per sequence:
##        seq No.mrk
##          1    411
##         10    245
##         11    284
##         12    221
##          2    281
##          3    318
##          4    432
##          5    291
##          6    433
##          7    409
##          8    279
##          9    310
##     ----------
##     Markers with no sequence information: 103
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  434
##      1  1  309
##      1  2  256
##      1  3  125
##      1  4   12
##      2  0  171
##      2  1  212
##      2  2  337
##      2  3  196
##      2  4   50
##      3  0   41
##      3  1  110
##      3  2  202
##      3  3  203
##      3  4  334
##      4  1    9
##      4  2   86
##      4  3  309
```

This function outputs information about the dataset including ploidy level, total number of individuals, total number of markers, number of informative markers and proportion of missing data. If `detailed = TRUE`, the function also outputs the number of markers in each sequence, if available, and the number of markers contained in all possible dosage combinations between both parents.

You can also explore the dataset visually using the `plot` function:


```r
plot(dat.dose.mpl)
```

![](README_files/figure-html/plot_data_dose-1.png)<!-- -->

The output figure shows a bar plot on the left-hand side with the number of markers contained in each allele dosage combination between both parents. The right labels indicate allele dosages for Parent 1 and Parent 2, respectively. The upper-right plot contains the $\log_10(p-value)$ from $\chi^2$ tests for to all markers, considering the expected segregation patterns under Mendelian inheritance. The lower-right plot contains a graphical representation of the allele dosages and missing data distribution for all markers and individuals.

### Marker-specific

If you want to view a specific marker information, use the `plot_mrk_info` function. You should indicate your dataset object using the `input.data` argument, and the desired marker using the `mrk` argument. You can indicate the marker using its number or its name (string):


```r
# For a dosage-based data analysis of marker 104
plot_mrk_info(input.data = dat.dose.mpl, mrk = 240)
```

![](README_files/figure-html/mrk_info-1.png)<!-- -->

```r
# For a probability-based data analysis of the marker solcap_snp_c1_13686
plot_mrk_info(input.data = dat.dist.mpl, mrk = 'solcap_snp_c1_13686')
```

![](README_files/figure-html/mrk_info-2.png)<!-- -->

When applied to a dosage-based dataset, the function outputs a figure showing: marker name and position in the dataset, allele dosage in parents 1 and 2, proportion of missing data, p-value of the associated $\chi^2$ test for Mendelian segregation, sequence and position information (when available). The figure also contains a plot with the allele dosage and missing data distibution in the population.

When applied to a probability-baed dataset, the function also outputs the probability threshold and a
3-dimensional plot containing the probability distribution for each allele dosage, considering all individuals.


## Filtering and Quality Control

In order to build a good genetic map, good quality data is desired to guarantee reliable estimates of recombination fractions and linkage phases. High proportions of messy data will reduce this reliability, such as unexpected segregation patterns and missing data. Furthermore, available technologies are able to produce hundreds of thousands of markers, and mid to high proportion of redundant markers is expected. Markers that carry the same information are expected to be positioned in the same place, and the redundant ones may be removed from dataset before starting the process, in order to reduce computational effort. MAPpoly handle some filtering functions, described in the sections below.

### Missing data filtering

MAPpoly is able to filter markers and/or individuals that exceeds a defined threshold for missing data. The function `filter_missing` does this and creates or updates your dataset following the rules passed by its arguments. The argument `input.data` should contain your dataset object, and you can choose to filter by 'marker' or 'individual' using the `type` argument. You can also define the maximum proportion of missing data using the `filter.thres` argument (range: 0-1, i.e. a threshold of 0.2 will keep just markers or individuals with less than 20% of missing data). When TRUE (default), the `inter` argument plots markers or individuals vs. frequency of missing data.
 

```r
# Filtering dataset by marker
dat.filt.mrk <- filter_missing(input.data = dat.dose.mpl, type = "marker", 
                               filter.thres = 0.2, inter = TRUE)
print(dat.filt.mrk)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             3978 
##     Missing data:                            2.7%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  432
##      1  1  301
##      1  2  254
##      1  3  123
##      1  4   12
##      2  0  170
##      2  1  205
##      2  2  333
##      2  3  196
##      2  4   49
##      3  0   41
##      3  1  107
##      3  2  200
##      3  3  199
##      3  4  333
##      4  1    9
##      4  2   86
##      4  3  307
```

```r
# Filtering dataset by individual
dat.filt.ind <- filter_missing(input.data = dat.filt.mrk, type = "individual", 
                               filter.thres = 0.1, inter = TRUE)
print(dat.filt.ind)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         154 
##     No. markers:                             3978 
##     Missing data:                            2.25%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  432
##      1  1  301
##      1  2  254
##      1  3  123
##      1  4   12
##      2  0  170
##      2  1  205
##      2  2  333
##      2  3  196
##      2  4   49
##      3  0   41
##      3  1  107
##      3  2  200
##      3  3  199
##      3  4  333
##      4  1    9
##      4  2   86
##      4  3  307
```

```r
dat.dose.filt <- dat.filt.ind
```



In this dataset, just 16 markers presented a proportion of missing data above the defined threshold, while no individual exceeded the defined threshold. Then we will keep the marker-based filtered dataset.

### Segregation test

Another very important point to consider is the expected marker segregation pattern under Mendelian inheritance. Markers with messy or distorted segregation would produce unreliable estimates, and may be removed (at least for a while) from the dataset. A good test for that is the chi-square ($\chi^2$) test, which basically matches expected genotype frequencies against observed frequencies and calculates the associated p-value. In order to define the p-value threshold for the tests, we will use the (conservative) Bonferroni correction: 

$$\alpha_{thres} = \frac{\alpha}{\#markers}$$

We will also assume that only random chromosome bivalent pairing occurs, and there is no double reduction.


```r
pval.bonf <- 0.05/dat.dose.filt$n.mrk
mrks.chi.filt <- filter_segregation(dat.dose.filt, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init<-make_seq_mappoly(mrks.chi.filt)
```

Please notice that `filter_segregation` does not produce a filtered dataset; it just tells you which markers follow the expected Mendelian segregation pattern. To select these markers from your dataset, you may use the `make_seq_mappoly` function.


```r
plot(seq.init)
```

![](README_files/figure-html/plot_filt_seq-1.png)<!-- -->


### Redundant markers

As mentioned before, mid to high proportion of redundant markers is expected from modern genotyping platforms. As markers that have the same information does not provide any advantage during the mapping process, redundant ones may be removed from the dataset in order to reduce computational effort. The function `elim_redundant` identifies which markers have the same information, and outputs the informative markers that should be used during the mapping process.


```r
seq.redundant <- elim_redundant(input.seq = seq.init)
print(seq.redundant)
```

```
## This is an object of class 'mappoly.sequence'
##     ------------------------
##     Parameters not estimated
##     ------------------------
##     Ploidy level:       4 
##     No. individuals:    154 
##     No. markers:        3598 
## 
##     ----------
##     No. markers per sequence:
##   sequence No.mrk
##          1    384
##         10    223
##         11    247
##         12    209
##          2    263
##          3    289
##          4    371
##          5    266
##          6    384
##          7    337
##          8    260
##          9    276
## 
##     ----------
##     No. of markers per dosage in both parents:
##     dP dQ freq
##      0  1  324
##      0  2  168
##      0  3   24
##      1  0  330
##      1  1  297
##      1  2  252
##      1  3  118
##      1  4   12
##      2  0  158
##      2  1  205
##      2  2  328
##      2  3  191
##      2  4   48
##      3  0   40
##      3  1  106
##      3  2  199
##      3  3  198
##      3  4  262
##      4  1    9
##      4  2   84
##      4  3  245
## ------------
## Eliminated markers:  379
```

Here you can see which markers are redundant and non-informative to the mapping process. Once you have identified them, you can represent the proportion graphically with:


```r
plot(seq.redundant)
```

![](README_files/figure-html/plot_redundant-1.png)<!-- -->

Please notice that `elim_redundant` does not produce a filtered sequence; it just tells you which markers are informative and which are redundant. To select just the informative ones, you may use the `make_seq_mappoly` function:


```r
seq.unique <- make_seq_mappoly(seq.redundant)
```

# Two-point analysis

Once the markers where selected, we need to compute the pairwise recombination fraction between all of them (two-point analysis). First, let us load the genotype counts ($\zeta_{\mbox{T}_{k},\mbox{T}_{k^{\prime}}}(l_{P}, l_{Q})$) defined in equation 20 in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). This object is fundamental to perform the dimension reduction of the transition space.


```r
counts<-cache_counts_twopt(input.seq = seq.unique, get.from.web = TRUE)
```

```
## Internet conectivety ok.
## Loading genotype counts from web
```

```r
counts
```

```
##   This is an object of class 'cache.info'
##   -----------------------------------------------------
##   Ploidy level:                                4 
##   No. marker combinations:                     625 
##   -----------------------------------------------------
```

The function `est_pairwise_rf` estimates all the pairwise recombination fractions in the sequence provided. Since the output object is too big to be fully displayed on the screen, `mappoly` shows a summary. Notice that parallel computation is available and in this case, we used 16 CPU's to perform the computations.



```r
#(~ 9.5 minutes)
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.unique, 
                                   count.cache = counts, 
                                   n.clusters = 16)
all.rf.pairwise
```

To assess the recombination fraction between a particular pair of markers, say 802 and 959, we use the following syntax:


```r
all.rf.pairwise$pairwise$`93-98`
```

```
##        LOD_ph         rf       LOD_rf
## 1-1  0.000000 0.06558141 3.4196648802
## 2-1 -3.083814 0.36564259 1.5102601666
## 1-0 -3.420131 0.49995416 0.0004662803
## 2-0 -3.540557 0.22902314 1.0535175904
## 0-1 -4.594369 0.49995416 0.0002950534
## 0-0 -4.595798 0.49992162 0.0017234837
```

```r
plot(all.rf.pairwise, first.mrk = 93, second.mrk = 98)
```

<img src="README_files/figure-html/plot_twopt_example-1.png" width="500px" style="display: block; margin: auto;" />

In this case, `93-98` represents the position of the markers in the filtered data set. The name of the rows have the form `x-y`, where `x` and `y` indicate how many homologous chromosomes share the same allelic variant in parents $P1$ and $P2$, respectively (see [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378) for notation). The first column indicates the LOD Score in relation to the most likely linkage phase configuration. The second column shows the recombination fraction, and the third indicates the LOD Score comparing the likelihood under no linkage ($r = 0.5$) and the estimated recombination fraction (evidence of linkage).

## Assembling the recombination fraction and LOD Score matrices

Recombination fraction and LOD Score matrices are fundamental in genetic mapping. Later in this tutorial, we will use these matrices as the basic information to order markers and also to perform some diagnostics. To convert the two-point object into recombination fraction and LOD Score matrices, we need to assume thresholds for the three columns observed in the previous output. The arguments `thresh.LOD.ph` and `thresh.LOD.rf` set LOD Scores thresholds for the second most likely linkage phase configuration and recombination fraction. Here we assume `thresh.LOD.ph = 0` and `thresh.LOD.rf = 0`, thus no matter how likely is the second best option, all the computed values will be considered. The argument `thresh.rf = 0.5` indicates that the maximum accepted recombination fraction is `0.5`. To convert these values in a recombination fraction matrix, we use the function `rf_list_to_matrix`


```r
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
```

```
## INFO: Going singlemode. Using one CPU.
```

We can plot this matrix using the reference genome order. For doing so, we use the function `get_genomic_order` to get the genomic order of the input sequence and use the resulting order to index the recombination fraction matrix. If the reference order is consistent with the marker order in this specific population, we should observe a block-diagonal matrix and within each sub-matrix, a monotonic pattern. 


```r
id<-get_genomic_order(seq.unique)
plot(mat, ord = rownames(id), index = FALSE)
```

<div class="figure" style="text-align: center">
<img src="image/all_rec_mat.png" alt="Figure 1: Example of CSV data set" width="90%" />
<p class="caption">Figure 1: Example of CSV data set</p>
</div>

As expected, we can observe the block-diagonal and monotonic patterns. In the previous case, the thresholds allowed to plot almost all points in the recombination fraction matrix. The empty cells in the matrix indicate markers where it is impossible to detect recombinant events using two-point estimates (e.g., between $1 \times 0$ and $0 \times 1$ marker). Yet, if the thresholds become more stringent (higher LODs and lower rf), the matrix becomes more sparse.  

# Assembling linkage groups

The function `group_mappoly` assign markers to linkage groups using the recombination fraction matrix obtained above. The user can provide an expected number of groups or run the interactive version of the function using `inter = TRUE`. Since in this data set we expect 12 linkage groups (basic chromosome number in potato), we use `expected.groups = 12`. If the data set provides the chromosome where the markers are located, the function allows comparing the groups obtained using the pairwise recombination fraction and the chromosome information provided using the `comp.mat = TRUE`.


```r
grs <- group_mappoly(input.mat = mat,
                     input.seq = seq.unique,
                     expected.groups = 12,
                     comp.mat = TRUE, 
                     inter = TRUE)
grs
```

```
##   This is an object of class 'mappoly.group'
##   ------------------------------------------
##   Criteria used to assign markers to groups:
## 
##     - Number of markers:          3598 
##     - Number of linkage groups:   12 
##     - Number of markers per linkage groups: 
##     group n.mrk
##         1   245
##         2   399
##         3   276
##         4   353
##         5   274
##         6   393
##         7   265
##         8   413
##         9   269
##        10   227
##        11   202
##        12   282
##   ------------------------------------------
##     11   6   9   7   5   4   2   1   8  10  12   3 NoChr
## 1  230   2   0   0   1   0   1   0   0   2   4   2     3
## 2    4 370   3   2   3   0   0   0   5   0   2   0    10
## 3    1   0 264   1   1   0   0   0   1   1   0   0     7
## 4    0   1   0 326   1   5   0   1   3   4   0   1    11
## 5    1   1   2   0 252   0   1   1   1   1   0   1    13
## 6    3   3   1   4   2 361   0   4   1   1   1   1    11
## 7    0   1   1   0   0   2 247   0   0   0   4   5     5
## 8    3   1   1   2   6   2  12 372   0   0   1   2    11
## 9    0   4   1   1   0   0   0   1 249   4   0   0     9
## 10   3   0   2   1   0   0   0   4   0 210   1   0     6
## 11   2   1   1   0   0   0   1   0   0   0 195   0     2
## 12   0   0   0   0   0   1   1   1   0   0   1 277     1
##   ------------------------------------------
```

Here, we have the 3598 markers distributed in 12 linkage groups. The rows indicate linkage groups obtained using linkage information and the columns are the chromosomes in the reference genome. Notice the diagonal indicating the concordance between the two sources of information. Now, we can plot the resulting marker cluster analysis.


```r
plot(grs)
```

![](README_files/figure-html/plot_group-1.png)<!-- -->

Once the linkage groups are properly assembled, we use the function `make_seq_mappoly` to make marker sequences from the group analysis. We will assemble a list with 12 positions, each one containing the corresponding linkage group sequence. Also, we will use only markers allocated in the diagonal of the previous comparison matrix. Thus only markers that were both assigned to a particular linkage group using both sources of information will be considered. We also will assemble smaller two-point objects to facilitate further parallelization procedures.


```r
LGS<-vector("list", 12)
for(j in 1:12){
  temp1<-make_seq_mappoly(grs, j)
  temp2<-get_genomic_order(temp1) # assembling sequence considering the genomic order
  lg.id<-as.numeric(names(which.max(table(temp2[,1]))))
  nm <- rownames(temp2)[which(temp2[,1] == lg.id)]
  temp3 <- make_seq_mappoly(dat.dose.filt, nm)
  tpt <- make_pairs_mappoly(all.rf.pairwise, input.seq = temp3)
  lgtemp<-rf_snp_filter(input.twopt = tpt)
  LGS[[lg.id]]<-list(lg = lgtemp, 
                 tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}
```

Now, let us print the recombination fraction matrices or each linkage group.

![](README_files/figure-html/all_mat_rf-1.png)<!-- -->

# Estimating the map for a given order

In this section, we will use the marker order provided by the _Solanum tuberosum_ genome version 4.03. The _de novo_ marker ordering will be addressed later in this tutorial. The estimation of the genetic map for a given order involves the computation of recombination fraction between adjacent markers and also finding the linkage phase configuration of those markers in both parents. The core function to perform these tasks in `mappoly` is `est_rf_hmm_sequential`. This function uses the pairwise recombination fraction as the first source of information to sequentially position allelic variants in specific homologs. For situations where pairwise analysis has limited power, the algorithm relies on the likelihood obtained through a hidden Markov model (HMM) [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378).  Once all markers are positioned, the final map is reconstructed using the HMM multipoint algorithm. 

![Example of linkage phase configuration estimation using sequential search space
reduction and HMM evaluation.](phasing.png)

Several arguments are available to control the inclusion and phasing of the markers in the chain.  `thres.twopt` receives the threshold to whether when the linkage phases compared via two-point analysis should be considered, and the HMM analysis should not be used to infer the linkage phase (A. K. A. $\eta$ in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378)). `thres.hmm` receives the threshold for keeping competing maps computed using HMM (if the two-point analysis was not enough) in the next round of marker insertion. `extend.tail` indicates the number of markers that should be considered at the end of the chain to insert a new marker. `tol` and `tol.final` receive the desired accuracy to estimate the sub-maps during the sequential phasing procedure and the desired accuracy in the final map. `phase.number.limit` receives the limit number of linkage phase configurations to be tested using HMM. `info.tail` is a logical argument and if `TRUE` it uses the complete informative tail (last markers in the chain that allow all homologous to be distinguished in the parents) of the chain to calculate the likelihood of the linkage phases. 

First, as an example, let us estimate the map for linkage group 3. The values used in the function arguments, were obtained using a balance of processing speed and accuracy of the algorithm. As an exercise, it is interesting to try different values and check out the results. For now, let us stick with the following values.


```r
lg3.map<-est_rf_hmm_sequential(input.seq = LGS[[3]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = LGS[[3]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```

```
## Number of markers: 263
```

```
## ═══════════════════════════════════════════════════════════════════════════════════════ Initial sequence ══
```

```
## 10 markers...
## ●    Trying sequence: 1 2 3 4 5 6 7 8 9 10 :
##        4 phase(s): . . . .
```

```
## ═════════════════════════════════════════════════════════════════════════════ Done with initial sequence ══
```

```
## 11: 4.2% (782): 2 ph(s) : (1/2)--t: 11
##     2.23: (-0.32/0.16)✔
## 12: 4.6% (783): 1 ph(s) : (1/1)--t: 12
##     3.1: (0.87/0.34)✔
## 13: 4.9% (786): 4 ph(s) : (1/4)--t: 13
##     3.99: (0.88/0.04)✔
## 14: 5.3% (787): 3 ph(s) : (1/3)--t: 14
##     10.09: (6.11/1.57)✔
## 15: 5.7% (788): 1 ph(s) : (1/1)--t: 15
##     10.09: (-0.01/-2.09)✔
## 16: 6.1% (789): 1 ph(s) : (1/1)--t: 16
##     10.3: (0.21/0.01)✔
## 17: 6.5% (790): 1 ph(s) : (1/1)--t: 17
##     10.55: (0.25/0.23)✔
## 18: 6.8% (791): 1 ph(s) : (1/1)--t: 18
##     10.56: (0.01/-0.24)✔
## 19: 7.2% (792): 1 ph(s) : (1/1)--t: 19
##     10.58: (0.01/0)✔
## 20: 7.6% (793): 1 ph(s) : (1/1)--t: 20
##     10.82: (0.24/0.08)✔
## 21: 8% (794): 1 ph(s) : (1/1)--t: 21
##     12.49: (1.66/1.57)✔
## 22: 8.4% (795): 1 ph(s) : (1/1)--t: 22
##     13.19: (0.7/-1.06)✔
## 23: 8.7% (796): 4 ph(s) : (1/4)--t: 23
##     17.34: (4.15/1.91)✔
## 24: 9.1% (797): 1 ph(s) : (1/1)--t: 24
##     18.35: (1.01/-1.47)✔
## 25: 9.5% (800): 1 ph(s) : (1/1)--t: 25
##     19.28: (0.93/-0.14)✔
## 26: 9.9% (801): 1 ph(s) : (1/1)--t: 26
##     20.81: (1.53/0.24)✔
## 27: 10.3% (802): 1 ph(s) : (1/1)--t: 27
##     21.56: (0.75/-0.67)✔
## 28: 10.6% (803): 16 ph(s) : (1/16)--t: 28
##     24.72: (3.16/0.87)
##     59.19: (37.62/30.76)
##     451.21: (429.65/425.38)
##     452: (430.44/425.38)✔
## 29: 11% (804): 4 ph(s) : (1/4)--t: 29
##     25.19: (0.47/-1.11)✔
## 30: 11.4% (805): 2 ph(s) : (1/2)--t: 30
##     26.11: (0.92/0.19)✔
## 31: 11.8% (806): 1 ph(s) : (1/1)--t: 31
##     26.09: (-0.02/-0.27)✔
## 32: 12.2% (807): 1 ph(s) : (1/1)--t: 32
##     26.14: (0.05/-0.11)✔
## 33: 12.5% (808): 1 ph(s) : (1/1)--t: 33
##     27.21: (1.07/0.44)✔
## 34: 12.9% (809): 2 ph(s) : (1/2)--t: 34
##     27.18: (-0.03/-0.47)✔
## 35: 13.3% (810): 1 ph(s) : (1/1)--t: 35
##     27.19: (0.01/0.03)✔
## 36: 13.7% (811): 1 ph(s) : (1/1)--t: 36
##     27.54: (0.35/0.12)✔
## 37: 14.1% (812): 1 ph(s) : (1/1)--t: 37
##     27.69: (0.15/-0.12)✔
## 38: 14.4% (813): 1 ph(s) : (1/1)--t: 38
##     27.71: (0.02/-0.04)✔
## 39: 14.8% (814): 1 ph(s) : (1/1)--t: 39
##     27.71: (0/0.03)✔
## 40: 15.2% (815): 1 ph(s) : (1/1)--t: 40
##     27.77: (0.06/0)✔
## 41: 15.6% (816): 1 ph(s) : (1/1)--t: 41
##     27.9: (0.13/-0.01)✔
## 42: 16% (818): 1 ph(s) : (1/1)--t: 42
##     28.45: (0.55/0.49)✔
## 43: 16.3% (819): 1 ph(s) : (1/1)--t: 43
##     29.09: (0.64/-0.06)✔
## 44: 16.7% (820): 1 ph(s) : (1/1)--t: 44
##     29.78: (0.69/0.11)✔
## 45: 17.1% (822): 1 ph(s) : (1/1)--t: 45
##     30.55: (0.77/-0.03)✔
## 46: 17.5% (823): 4 ph(s) : (1/4)--t: 46
##     30.39: (-0.16/-0.5)✔
## 47: 17.9% (824): 1 ph(s) : (1/1)--t: 47
##     31.33: (0.94/0.45)✔
## 48: 18.3% (825): 1 ph(s) : (1/1)--t: 48
##     31.27: (-0.06/-0.45)✔
## 49: 18.6% (826): 1 ph(s) : (1/1)--t: 49
##     32.6: (1.33/0.53)✔
## 50: 19% (827): 1 ph(s) : (1/1)--t: 50
##     32.57: (-0.04/-0.53)✔
## 51: 19.4% (828): 1 ph(s) : (1/1)--t: 51
##     32.77: (0.2/0.1)✔
## 52: 19.8% (829): 1 ph(s) : (1/1)--t: 52
##     32.82: (0.06/-0.1)✔
## 53: 20.2% (836): 1 ph(s) : (1/1)--t: 53
##     32.85: (0.03/0.04)✔
## 54: 20.5% (837): 1 ph(s) : (1/1)--t: 54
##     33.14: (0.28/0.11)✔
## 55: 20.9% (838): 1 ph(s) : (1/1)--t: 55
##     34.17: (1.04/0.29)✔
## 56: 21.3% (840): 1 ph(s) : (1/1)--t: 56
##     34.22: (0.04/-0.43)✔
## 57: 21.7% (841): 1 ph(s) : (1/1)--t: 57
##     34.45: (0.23/0.01)✔
## 58: 22.1% (846): 1 ph(s) : (1/1)--t: 58
##     34.5: (0.06/-0.01)✔
## 59: 22.4% (847): 1 ph(s) : (1/1)--t: 59
##     34.47: (-0.04/-0.03)✔
## 60: 22.8% (848): 1 ph(s) : (1/1)--t: 60
##     34.79: (0.33/0.23)✔
## 61: 23.2% (849): 4 ph(s) : (1/4)--t: 61
##     35.19: (0.4/-0.09)✔
## 62: 23.6% (850): 1 ph(s) : (1/1)--t: 62
##     35.98: (0.79/0.62)✔
## 63: 24% (851): 1 ph(s) : (1/1)--t: 63
##     37.13: (1.16/0.07)✔
## 64: 24.3% (852): 1 ph(s) : (1/1)--t: 64
##     37.65: (0.52/-0.34)✔
## 65: 24.7% (853): 1 ph(s) : (1/1)--t: 65
##     37.65: (0/-0.5)✔
## 66: 25.1% (854): 1 ph(s) : (1/1)--t: 66
##     37.92: (0.27/0.24)✔
## 67: 25.5% (855): 3 ph(s) : (1/3)--t: 67
##     38.45: (0.53/0.06)✔
## 68: 25.9% (856): 1 ph(s) : (1/1)--t: 68
##     38.69: (0.24/-0.08)✔
## 69: 26.2% (857): 1 ph(s) : (1/1)--t: 69
##     38.71: (0.02/-0.23)✔
## 70: 26.6% (858): 1 ph(s) : (1/1)--t: 70
##     38.95: (0.25/0.06)✔
## 71: 27% (859): 1 ph(s) : (1/1)--t: 71
##     38.77: (-0.18/0.02)✔
## 72: 27.4% (860): 1 ph(s) : (1/1)--t: 72
##     38.83: (0.06/-0.04)✔
## 73: 27.8% (861): 1 ph(s) : (1/1)--t: 73
##     39.15: (0.32/0.21)✔
## 74: 28.1% (862): 1 ph(s) : (1/1)--t: 74
##     39.38: (0.23/-0.08)✔
## 75: 28.5% (863): 1 ph(s) : (1/1)--t: 75
##     40.32: (0.94/0.08)✔
## 76: 28.9% (864): 1 ph(s) : (1/1)--t: 76
##     40.42: (0.1/-0.09)✔
## 77: 29.3% (865): 1 ph(s) : (1/1)--t: 77
##     41.97: (1.55/1.21)✔
## 78: 29.7% (866): 1 ph(s) : (1/1)--t: 78
##     42.38: (0.41/-1.11)✔
## 79: 30% (867): 1 ph(s) : (1/1)--t: 79
##     42.47: (0.09/-0.07)✔
## 80: 30.4% (868): 1 ph(s) : (1/1)--t: 80
##     43.41: (0.94/0.24)✔
## 81: 30.8% (869): 1 ph(s) : (1/1)--t: 81
##     44.02: (0.61/0.08)✔
## 82: 31.2% (870): 1 ph(s) : (1/1)--t: 82
##     44.24: (0.22/-0.44)✔
## 83: 31.6% (871): 1 ph(s) : (1/1)--t: 83
##     44.5: (0.26/0.19)✔
## 84: 31.9% (872): 1 ph(s) : (1/1)--t: 84
##     44.9: (0.4/-0.09)✔
## 85: 32.3% (873): 1 ph(s) : (1/1)--t: 85
##     45.04: (0.14/-0.07)✔
## 86: 32.7% (874): 1 ph(s) : (1/1)--t: 86
##     46.1: (1.06/0.93)✔
## 87: 33.1% (875): 1 ph(s) : (1/1)--t: 87
##     46.25: (0.15/-0.85)✔
## 88: 33.5% (876): 1 ph(s) : (1/1)--t: 88
##     46.91: (0.66/0.17)✔
## 89: 33.8% (877): 1 ph(s) : (1/1)--t: 89
##     47.4: (0.48/-0.15)✔
## 90: 34.2% (878): 1 ph(s) : (1/1)--t: 90
##     48.26: (0.86/0.31)✔
## 91: 34.6% (879): 1 ph(s) : (1/1)--t: 91
##     48.73: (0.47/-0.02)✔
## 92: 35% (880): 1 ph(s) : (1/1)--t: 92
##     48.7: (-0.03/-0.48)✔
## 93: 35.4% (881): 1 ph(s) : (1/1)--t: 93
##     49.83: (1.13/0.35)✔
## 94: 35.7% (883): 1 ph(s) : (1/1)--t: 94
##     49.82: (-0.01/-0.35)✔
## 95: 36.1% (884): 1 ph(s) : (1/1)--t: 95
##     50.77: (0.95/0.9)✔
## 96: 36.5% (885): 1 ph(s) : (1/1)--t: 96
##     51.13: (0.36/-0.76)✔
## 97: 36.9% (886): 1 ph(s) : (1/1)--t: 97
##     51.39: (0.27/0.1)✔
## 98: 37.3% (887): 1 ph(s) : (1/1)--t: 98
##     52.16: (0.77/0.5)✔
## 99: 37.6% (888): 1 ph(s) : (1/1)--t: 99
##     52.47: (0.31/-0.64)✔
## 100: 38% (889): 1 ph(s) : (1/1)--t: 100
##     52.79: (0.32/0.05)✔
## 101: 38.4% (890): 1 ph(s) : (1/1)--t: 101
##     53.16: (0.37/0.01)✔
## 102: 38.8% (891): 1 ph(s) : (1/1)--t: 102
##     53.62: (0.46/0.14)✔
## 103: 39.2% (892): 1 ph(s) : (1/1)--t: 103
##     54.38: (0.77/0.24)✔
## 104: 39.5% (893): 1 ph(s) : (1/1)--t: 104
##     55.82: (1.44/0.95)✔
## 105: 39.9% (894): 1 ph(s) : (1/1)--t: 105
##     56.04: (0.22/-1.28)✔
## 106: 40.3% (896): 1 ph(s) : (1/1)--t: 106
##     57.63: (1.58/0.72)✔
## 107: 40.7% (898): 1 ph(s) : (1/1)--t: 107
##     57.75: (0.12/-0.82)✔
## 108: 41.1% (899): 1 ph(s) : (1/1)--t: 108
##     57.62: (-0.12/-0.08)✔
## 109: 41.4% (900): 1 ph(s) : (1/1)--t: 109
##     58.14: (0.52/0.53)✔
## 110: 41.8% (902): 1 ph(s) : (1/1)--t: 110
##     58.18: (0.04/-0.55)✔
## 111: 42.2% (903): 1 ph(s) : (1/1)--t: 111
##     58.01: (-0.17/0.01)✔
## 112: 42.6% (906): 1 ph(s) : (1/1)--t: 112
##     58.95: (0.95/0.47)✔
## 113: 43% (907): 1 ph(s) : (1/1)--t: 113
##     59.2: (0.24/-0.36)✔
## 114: 43.3% (908): 1 ph(s) : (1/1)--t: 114
##     59.3: (0.1/-0.1)✔
## 115: 43.7% (909): 1 ph(s) : (1/1)--t: 115
##     60.98: (1.69/1.06)✔
## 116: 44.1% (910): 1 ph(s) : (1/1)--t: 116
##     62.25: (1.26/0.02)✔
## 117: 44.5% (911): 1 ph(s) : (1/1)--t: 117
##     63.09: (0.85/-0.23)✔
## 118: 44.9% (912): 1 ph(s) : (1/1)--t: 118
##     63.15: (0.06/-0.86)✔
## 119: 45.2% (914): 1 ph(s) : (1/1)--t: 119
##     64.1: (0.95/0.25)✔
## 120: 45.6% (915): 1 ph(s) : (1/1)--t: 120
##     64.05: (-0.04/-0.21)✔
## 121: 46% (916): 1 ph(s) : (1/1)--t: 121
##     65.15: (1.1/0.67)✔
## 122: 46.4% (917): 1 ph(s) : (1/1)--t: 122
##     65.32: (0.17/-0.63)✔
## 123: 46.8% (918): 1 ph(s) : (1/1)--t: 123
##     65.57: (0.25/0.05)✔
## 124: 47.1% (919): 2 ph(s) : (1/2)--t: 124
##     66.25: (0.68/0.18)✔
## 125: 47.5% (920): 1 ph(s) : (1/1)--t: 125
##     66.57: (0.32/-0.06)✔
## 126: 47.9% (921): 1 ph(s) : (1/1)--t: 126
##     66.8: (0.24/-0.21)✔
## 127: 48.3% (922): 1 ph(s) : (1/1)--t: 127
##     67.65: (0.85/0.43)✔
## 128: 48.7% (924): 1 ph(s) : (1/1)--t: 128
##     67.9: (0.25/-0.21)✔
## 129: 49% (925): 1 ph(s) : (1/1)--t: 129
##     67.77: (-0.13/-0.26)✔
## 130: 49.4% (926): 1 ph(s) : (1/1)--t: 130
##     67.83: (0.06/0.09)✔
## 131: 49.8% (927): 1 ph(s) : (1/1)--t: 131
##     68.49: (0.66/0.43)✔
## 132: 50.2% (928): 1 ph(s) : (1/1)--t: 132
##     68.65: (0.16/-0.24)✔
## 133: 50.6% (929): 1 ph(s) : (1/1)--t: 133
##     68.84: (0.18/-0.27)✔
## 134: 51% (930): 1 ph(s) : (1/1)--t: 134
##     69.42: (0.58/0.57)✔
## 135: 51.3% (931): 1 ph(s) : (1/1)--t: 135
##     70.57: (1.16/-0.3)✔
## 136: 51.7% (932): 1 ph(s) : (1/1)--t: 136
##     71.11: (0.53/-0.03)✔
## 137: 52.1% (933): 1 ph(s) : (1/1)--t: 137
##     72.43: (1.32/0.41)✔
## 138: 52.5% (934): 1 ph(s) : (1/1)--t: 138
##     72.42: (-0.02/-0.66)✔
## 139: 52.9% (935): 1 ph(s) : (1/1)--t: 139
##     73.53: (1.11/0.45)✔
## 140: 53.2% (936): 1 ph(s) : (1/1)--t: 140
##     74.01: (0.48/-0.19)✔
## 141: 53.6% (937): 1 ph(s) : (1/1)--t: 141
##     74.43: (0.42/-0.14)✔
## 142: 54% (938): 1 ph(s) : (1/1)--t: 142
##     74.67: (0.24/-0.06)✔
## 143: 54.4% (939): 1 ph(s) : (1/1)--t: 143
##     74.78: (0.11/0.16)✔
## 144: 54.8% (940): 1 ph(s) : (1/1)--t: 144
##     75.12: (0.35/0)✔
## 145: 55.1% (941): 1 ph(s) : (1/1)--t: 145
##     75.14: (0.02/-0.25)✔
## 146: 55.5% (942): 1 ph(s) : (1/1)--t: 146
##     75.16: (0.01/0)✔
## 147: 55.9% (943): 1 ph(s) : (1/1)--t: 147
##     75.18: (0.02/0)✔
## 148: 56.3% (944): 1 ph(s) : (1/1)--t: 148
##     75.47: (0.29/0.32)✔
## 149: 56.7% (945): 1 ph(s) : (1/1)--t: 149
##     75.62: (0.15/-0.08)✔
## 150: 57% (946): 1 ph(s) : (1/1)--t: 150
##     76.74: (1.12/0.17)✔
## 151: 57.4% (947): 1 ph(s) : (1/1)--t: 151
##     77.2: (0.46/0.04)✔
## 152: 57.8% (948): 1 ph(s) : (1/1)--t: 152
##     77.11: (-0.09/-0.45)✔
## 153: 58.2% (949): 1 ph(s) : (1/1)--t: 153
##     77.1: (0/-0.01)✔
## 154: 58.6% (950): 1 ph(s) : (1/1)--t: 154
##     77.34: (0.24/0.23)✔
## 155: 58.9% (951): 1 ph(s) : (1/1)--t: 155
##     77.58: (0.23/0.1)✔
## 156: 59.3% (952): 1 ph(s) : (1/1)--t: 156
##     77.95: (0.38/-0.11)✔
## 157: 59.7% (953): 1 ph(s) : (1/1)--t: 157
##     77.34: (-0.61/0.13)✔
## 158: 60.1% (954): 1 ph(s) : (1/1)--t: 158
##     78.42: (1.08/0.23)✔
## 159: 60.5% (955): 2 ph(s) : (1/2)--t: 159
##     78.97: (0.54/-0.12)✔
## 160: 60.8% (956): 1 ph(s) : (1/1)--t: 160
##     79.53: (0.57/-0.04)✔
## 161: 61.2% (957): 1 ph(s) : (1/1)--t: 161
##     82.22: (2.69/2.26)✔
## 162: 61.6% (958): 1 ph(s) : (1/1)--t: 162
##     82.3: (0.08/-2.35)✔
## 163: 62% (959): 1 ph(s) : (1/1)--t: 163
##     83.02: (0.72/0.43)✔
## 164: 62.4% (960): 1 ph(s) : (1/1)--t: 164
##     84.23: (1.21/-0.49)✔
## 165: 62.7% (961): 1 ph(s) : (1/1)--t: 165
##     86.23: (2/0.25)✔
## 166: 63.1% (962): 1 ph(s) : (1/1)--t: 166
##     86.74: (0.51/-0.12)✔
## 167: 63.5% (963): 1 ph(s) : (1/1)--t: 167
##     86.77: (0.04/-0.36)✔
## 168: 63.9% (964): 1 ph(s) : (1/1)--t: 168
##     88.03: (1.25/0.36)✔
## 169: 64.3% (965): 1 ph(s) : (1/1)--t: 169
##     89.27: (1.24/0.71)✔
## 170: 64.6% (968): 1 ph(s) : (1/1)--t: 170
##     90.28: (1.02/-0.42)✔
## 171: 65% (969): 1 ph(s) : (1/1)--t: 171
##     90.82: (0.54/-0.33)✔
## 172: 65.4% (970): 1 ph(s) : (1/1)--t: 172
##     90.86: (0.03/-0.22)✔
## 173: 65.8% (973): 1 ph(s) : (1/1)--t: 173
##     91.12: (0.26/0.03)✔
## 174: 66.2% (974): 3 ph(s) : (1/3)--t: 174
##     92.23: (1.12/0.3)✔
## 175: 66.5% (975): 1 ph(s) : (1/1)--t: 175
##     92.49: (0.25/-0.16)✔
## 176: 66.9% (976): 1 ph(s) : (1/1)--t: 176
##     92.9: (0.41/0.2)✔
## 177: 67.3% (977): 1 ph(s) : (1/1)--t: 177
##     92.9: (0/-0.5)✔
## 178: 67.7% (978): 1 ph(s) : (1/1)--t: 178
##     93.57: (0.67/0.29)✔
## 179: 68.1% (979): 1 ph(s) : (1/1)--t: 179
##     94.14: (0.57/0.02)✔
## 180: 68.4% (980): 1 ph(s) : (1/1)--t: 180
##     95.18: (1.04/0.47)✔
## 181: 68.8% (983): 1 ph(s) : (1/1)--t: 181
##     97.16: (1.99/0.62)✔
## 182: 69.2% (984): 2 ph(s) : (1/2)--t: 182
##     120.38: (23.22/13.79)
##     121.73: (24.57/14.93)✖
## 182: not included (map expansion)
## 183: 69.6% (985): 1 ph(s) : (1/1)--t: 182
##     98.51: (1.35/-0.43)✔
## 184: 70% (986): 2 ph(s) : (1/2)--t: 183
##     121.59: (23.08/17.5)
##     125.04: (26.53/20.81)✖
## 184: not included (map expansion)
## 185: 70.3% (987): 1 ph(s) : (1/1)--t: 183
##     98.63: (0.11/-0.87)✔
## 186: 70.7% (988): 1 ph(s) : (1/1)--t: 184
##     98.66: (0.04/0.01)✔
## 187: 71.1% (989): 1 ph(s) : (1/1)--t: 185
##     100.63: (1.96/1.16)✔
## 188: 71.5% (990): 1 ph(s) : (1/1)--t: 186
##     101.16: (0.53/-0.65)✔
## 189: 71.9% (991): 1 ph(s) : (1/1)--t: 187
##     101.52: (0.36/-0.45)✔
## 190: 72.2% (992): 1 ph(s) : (1/1)--t: 188
##     102.72: (1.2/0.21)✔
## 191: 72.6% (993): 1 ph(s) : (1/1)--t: 189
##     103.09: (0.37/-0.11)✔
## 192: 73% (994): 1 ph(s) : (1/1)--t: 190
##     103.15: (0.06/-0.23)✔
## 193: 73.4% (996): 1 ph(s) : (1/1)--t: 191
##     103.32: (0.17/0.11)✔
## 194: 73.8% (997): 1 ph(s) : (1/1)--t: 192
##     104.08: (0.76/0.58)✔
## 195: 74.1% (998): 1 ph(s) : (1/1)--t: 193
##     104.29: (0.21/-0.61)✔
## 196: 74.5% (1000): 1 ph(s) : (1/1)--t: 194
##     105.93: (1.64/0.59)✔
## 197: 74.9% (1001): 1 ph(s) : (1/1)--t: 195
##     106.65: (0.71/-0.41)✔
## 198: 75.3% (1002): 1 ph(s) : (1/1)--t: 196
##     108.36: (1.71/0.74)✔
## 199: 75.7% (1003): 1 ph(s) : (1/1)--t: 197
##     109.62: (1.26/0.31)✔
## 200: 76% (1004): 1 ph(s) : (1/1)--t: 198
##     109.79: (0.17/-1.02)✔
## 201: 76.4% (1005): 1 ph(s) : (1/1)--t: 199
##     110.42: (0.63/0.28)✔
## 202: 76.8% (1006): 1 ph(s) : (1/1)--t: 200
##     110.47: (0.06/-0.57)✔
## 203: 77.2% (1007): 1 ph(s) : (1/1)--t: 201
##     110.82: (0.34/0.14)✔
## 204: 77.6% (1009): 1 ph(s) : (1/1)--t: 201
##     110.62: (-0.2/-0.14)✔
## 205: 77.9% (1012): 1 ph(s) : (1/1)--t: 201
##     109.54: (-1.08/0.27)✔
## 206: 78.3% (1013): 1 ph(s) : (1/1)--t: 201
##     110.4: (0.87/0.47)✔
## 207: 78.7% (1014): 1 ph(s) : (1/1)--t: 201
##     109.86: (-0.54/-0.74)✔
## 208: 79.1% (1015): 1 ph(s) : (1/1)--t: 201
##     110.3: (0.44/0.37)✔
## 209: 79.5% (1016): 2 ph(s) : (1/2)--t: 201
##     137.61: (27.31/24.56)
##     146.88: (36.58/33.99)✖
## 209: not included (map expansion)
## 210: 79.8% (1017): 1 ph(s) : (1/1)--t: 201
##     110.84: (0.54/-0.08)✔
## 211: 80.2% (1018): 1 ph(s) : (1/1)--t: 201
##     110.69: (-0.15/0)✔
## 212: 80.6% (1020): 1 ph(s) : (1/1)--t: 201
##     110.11: (-0.58/0.05)✔
## 213: 81% (1021): 1 ph(s) : (1/1)--t: 201
##     107.56: (-2.54/-0.25)✔
## 214: 81.4% (1023): 1 ph(s) : (1/1)--t: 201
##     109.35: (1.78/1.02)✔
## 215: 81.7% (1024): 1 ph(s) : (1/1)--t: 201
##     108.54: (-0.81/-0.65)✔
## 216: 82.1% (1025): 1 ph(s) : (1/1)--t: 201
##     109.1: (0.56/0.06)✔
## 217: 82.5% (1026): 1 ph(s) : (1/1)--t: 201
##     107.61: (-1.49/-0.42)✔
## 218: 82.9% (1027): 1 ph(s) : (1/1)--t: 201
##     107.48: (-0.13/-0.1)✔
## 219: 83.3% (1028): 1 ph(s) : (1/1)--t: 201
##     107.8: (0.33/0.28)✔
## 220: 83.7% (1029): 1 ph(s) : (1/1)--t: 201
##     109.33: (1.53/0.36)✔
## 221: 84% (1030): 1 ph(s) : (1/1)--t: 201
##     110.02: (0.69/-0.25)✔
## 222: 84.4% (1031): 1 ph(s) : (1/1)--t: 201
##     110.25: (0.23/-0.19)✔
## 223: 84.8% (1032): 1 ph(s) : (1/1)--t: 201
##     110.73: (0.48/0.23)✔
## 224: 85.2% (1033): 1 ph(s) : (1/1)--t: 201
##     109.11: (-1.61/-0.47)✔
## 225: 85.6% (1034): 1 ph(s) : (1/1)--t: 201
##     108.03: (-1.08/0.11)✔
## 226: 85.9% (1035): 2 ph(s) : (1/2)--t: 201
##     106.34: (-1.69/0.26)✔
## 227: 86.3% (1036): 1 ph(s) : (1/1)--t: 201
##     105.38: (-0.96/0.3)✔
## 228: 86.7% (1037): 1 ph(s) : (1/1)--t: 201
##     104.29: (-1.09/-0.52)✔
## 229: 87.1% (1038): 1 ph(s) : (1/1)--t: 201
##     103.86: (-0.43/-0.05)✔
## 230: 87.5% (1039): 1 ph(s) : (1/1)--t: 201
##     101.78: (-2.08/0.2)✔
## 231: 87.8% (1040): 1 ph(s) : (1/1)--t: 201
##     102.12: (0.34/0.5)✔
## 232: 88.2% (1041): 1 ph(s) : (1/1)--t: 201
##     102.79: (0.67/0.79)✔
## 233: 88.6% (1042): 1 ph(s) : (1/1)--t: 201
##     103.27: (0.48/-1.32)✔
## 234: 89% (1043): 1 ph(s) : (1/1)--t: 201
##     102.87: (-0.4/-0.21)✔
## 235: 89.4% (1045): 1 ph(s) : (1/1)--t: 201
##     102.89: (0.02/0.17)✔
## 236: 89.7% (1046): 1 ph(s) : (1/1)--t: 201
##     104.93: (2.03/1.83)✔
## 237: 90.1% (1047): 1 ph(s) : (1/1)--t: 201
##     105.2: (0.28/-1.85)✔
## 238: 90.5% (1048): 1 ph(s) : (1/1)--t: 201
##     105.11: (-0.1/-0.18)✔
## 239: 90.9% (1049): 1 ph(s) : (1/1)--t: 201
##     105.03: (-0.08/0)✔
## 240: 91.3% (1050): 1 ph(s) : (1/1)--t: 201
##     104.62: (-0.41/0)✔
## 241: 91.6% (1051): 1 ph(s) : (1/1)--t: 201
##     106.37: (1.75/0.43)✔
## 242: 92% (1052): 1 ph(s) : (1/1)--t: 201
##     108.24: (1.86/1.02)✔
## 243: 92.4% (1053): 1 ph(s) : (1/1)--t: 201
##     109.25: (1.01/-0.24)✔
## 244: 92.8% (1054): 1 ph(s) : (1/1)--t: 201
##     109.33: (0.08/-0.98)✔
## 245: 93.2% (1055): 3 ph(s) : (1/3)--t: 201
##     133.27: (23.95/20.93)
##     137.17: (27.84/24.71)
##     134.92: (25.6/22.58)✖
## 245: not included (map expansion)
## 246: 93.5% (1056): 1 ph(s) : (1/1)--t: 201
##     110.68: (1.35/0.36)✔
## 247: 93.9% (1057): 1 ph(s) : (1/1)--t: 201
##     110.37: (-0.3/-0.38)✔
## 248: 94.3% (1058): 1 ph(s) : (1/1)--t: 201
##     110.7: (0.33/0.18)✔
## 249: 94.7% (1059): 1 ph(s) : (1/1)--t: 201
##     109.99: (-0.71/-0.17)✔
## 250: 95.1% (1060): 1 ph(s) : (1/1)--t: 201
##     109.81: (-0.18/-0.05)✔
## 251: 95.4% (1063): 1 ph(s) : (1/1)--t: 201
##     108.63: (-1.19/-0.1)✔
## 252: 95.8% (1065): 1 ph(s) : (1/1)--t: 201
##     108.87: (0.25/0.18)✔
## 253: 96.2% (1066): 1 ph(s) : (1/1)--t: 201
##     109.25: (0.38/-0.01)✔
## 254: 96.6% (1067): 1 ph(s) : (1/1)--t: 201
##     109.17: (-0.08/-0.21)✔
## 255: 97% (1068): 1 ph(s) : (1/1)--t: 201
##     108.97: (-0.2/0.52)✔
## 256: 97.3% (1069): 1 ph(s) : (1/1)--t: 201
##     114.47: (5.5/2.69)✔
## 257: 97.7% (1070): 2 ph(s) : (1/2)--t: 201
##     115.73: (1.26/-1.47)
##     540.97: (426.5/422.55)✔
## 258: 98.1% (1071): 1 ph(s) : (1/1)--t: 201
##     117.45: (1.72/-1.79)✔
## 259: 98.5% (1073): 1 ph(s) : (1/1)--t: 201
##     123.18: (5.72/2.93)✔
## 260: 98.9% (1074): 1 ph(s) : (1/1)--t: 201
##     124.06: (0.89/-1.89)✔
## 261: 99.2% (1076): 1 ph(s) : (1/1)--t: 201
##     125.49: (1.43/-0.26)✔
## 262: 99.6% (1077): 1 ph(s) : (1/1)--t: 201
##     127.62: (2.13/-0.03)✔
## 263: 100% (1079): 1 ph(s) : (1/1)--t: 201
##     127.66: (0.04/-0.37)✔
```

```
## ═════════════════════════════════════════════════════════════ Reestimating final recombination fractions ══
```

```
## Markers in the initial sequence: 263
## Maped markers                  : 259 (98.5%)
```

```
## ═══════════════════════════════════════════════════════════════════════════════════════════════════════════
```

Now, we can display the results using the functions `print` and `plot`.


```r
print(lg3.map)
```

```
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 154 
##     No. markers:	 259 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 165.85
##        log-likelihood:	 -2271.88
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
```

```r
plot(lg3.map)
```

![](README_files/figure-html/map_lg3_plot-1.png)<!-- -->

Colored rectangles (red and blue) indicates the presence of the allelic variant in each one of the four homologous in both parents, $P_1$ and $P_2$.

## Reestimaing the map considering genotyping errors

Though current technologies enabled the genotyping of thousands of SNPs, they are quite prone to genotyping errors. One way to address this problem is to associate a probability distribution to each one of the markers and allow the HMM to update their probability. We can apply this procedure using either the probability distribution provided by the genotype calling software (loaded in MAPpoly using the function 'read_geno_dist') or assuming a global genotype error. For a detailed explanation of this procedure, please see [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378).  In this case, since we are analyzing the dosage data, with no probability distribution associated, we are going to use the second approch. Briefly, the use of the prior information will update the genotype of the markers based on a global chromosome structure. In this tutorial, since we are analyzing the dosage data, with no probability distribution associated, we are going to use the second approch. 


```r
lg3.map.error<-est_full_hmm_with_global_error(input.map = lg3.map, error = 0.05)
```

```
## 
##  ----------------------------------------------
##  INFO: running HMM using full transition space:
##        this operation may take a while.
## -----------------------------------------------
```

```r
plot(lg3.map.error)
```

![](README_files/figure-html/map_prior_lg3-1.png)<!-- -->

Notice that we considered a global genotyping error of 5% and the resulting map is smaller than the previous. Also, some markers were "attracted" and some markers were "repealed". 

# Posterior genotype probabilities

Now, let us compute the posterior genotype probabilities which is a peace of information necessary to perform further QTL analysis. Let us use the function `calc_genoprob_error`, whcih, similarly to `est_full_hmm_with_global_error` allows the inclusion of a global genotyping error 


```r
genoprob.lg3<-calc_genoprob_error(input.map = lg3.map, error = 0.05)
```

```
## Ploidy level:4
## Number of individuals:154
## 	..................................................
## 	..................................................
## 	..................................................
## 	....
```

Again, we used a lobal genotyping error of 5%. Each position of the object `genoprob.lg3` contains two elements: an array of dimensions $36 \times number \; of \; markers \times  number \; of \; individuals$ and the position of the markers in the maps in centimorgans.  Let us display the results for individual 1 in the full-sib popualtion


```r
ind <- 1
d <- genoprob.lg3$map
pr <- genoprob.lg3$probs[,,ind]
image(t(pr),
      col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
      axes=FALSE,
      xlab = "Markers",
      ylab = "",
      main = paste("LG", i))
axis(side = 1, at = d/max(d),
     labels =rep("", length(d)), las=2)
axis(side = 2, at = seq(0,1,length.out = nrow(pr)),
     labels = rownames(pr), las=2, cex.axis=.5)
```

![](README_files/figure-html/plot_genoprob_prior-1.png)<!-- -->

In this figure, the x-axis represents the genetic map and the y-axis represents the 36 possible genotypes in the full-sib population. The color scale varies from dark purple (high probabilityes) to light yellow (low probabilities). The `genoprob` object obtained here can be used to perform QTL analysis using the R package [QTLpoly](https://github.com/guilherme-pereira/QTLpoly), which is an under development software to map multiple QTLs in full-sib families of outcrossing autopolyploid species. 

# Ordering markers using MDS and reestimating the map

So far we reestimated the map using the genomic order. In real situations, unless a genomic information is provided, we need to order the markers using an optimization technique. Here, we use the MDS (multidimensional scaling) algorithm, proposed in the context of genetic mapping by [Preedy and Hackett (2016)](https://link.springer.com/article/10.1007%2Fs00122-016-2761-8). The MDS algorithm requires a recombination fraction matrix, which will be transformed in distance using a mapping function (in this case we use Haldane's mapping function). First, let us gather the pairwise recombination fractions for all three linkage groups


```r
mt <- lapply(LGS, function(x) rf_list_to_matrix(x$tpt))
```

Now, for each matrix contained in the object in `mt`, we use the MDS algorithm


```r
mds.ord <- lapply(mt, mds_mappoly)
```

```
## Registered S3 methods overwritten by 'MDSMap':
##   method       from   
##   plot.pcmap   mappoly
##   plot.pcmap3d mappoly
```
Usually, at this point, the user plot diagnostic plots to remove markers that are disturbing the ordering procedure. Here we didn't use that procedure, but we encourage the user to check the example in `?mds_mappoly`. Now, let us compare the estimated and the simulated order 

```r
LGS.mds<-vector("list", 12)
for(j in 1:12){
  lgtemp<-make_seq_mappoly(mds.ord[[j]])
  LGS.mds[[j]]<-list(lg = lgtemp, 
                 tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}
op <- par(mfrow = c(3, 4), pty = "s", mar=c(1,1,1,1)) 
sapply(LGS.mds, function(x) {
  plot(x = order(x$lg$sequence.pos), 
       y = rank(x$lg$seq.num), 
       xlab = "genomic order", 
       ylab = "estimated order")
  })
```

![](README_files/figure-html/compare_order-1.png)<!-- -->

```r
par(op)
```

Although we can observe several local inconsistencies, the global diagonal patterns indicate a very good order for all linkage groups. Now, let us build the genetic map of linkage group 3 using the MDS order


```r
lg3.map.mds<-est_rf_hmm_sequential(input.seq = LGS.mds[[3]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = LGS.mds[[3]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```

```
## Number of markers: 263
```

```
## ═══════════════════════════════════════════════════════════════════════════════════════ Initial sequence ══
```

```
## 10 markers...
## ●    Trying sequence: 1 2 3 4 5 6 7 8 9 10 :
##        28 phase(s): . . . . . . . . . . . . . . . . . . . . . . . . . . . .
```

```
## ═════════════════════════════════════════════════════════════════════════════ Done with initial sequence ══
```

```
## 11: 4.2% (1025): 2 ph(s) : (1/2)--t: 11
##     21.76: (1.28/-0.68)
##     446.43: (425.96/422.05)✔
## 12: 4.6% (1065): 1 ph(s) : (1/1)--t: 12
##     25.48: (3.72/-1.41)✔
## 13: 4.9% (1052): 1 ph(s) : (1/1)--t: 13
##     26.65: (1.18/-0.78)✔
## 14: 5.3% (1058): 1 ph(s) : (1/1)--t: 14
##     26.67: (0.02/-0.92)✔
## 15: 5.7% (1059): 4 ph(s) : (1/4)--t: 15
##     26.13: (-0.54/0.11)✔
## 16: 6.1% (1034): 3 ph(s) : (1/3)--t: 16
##     28.54: (2.41/2.26)
##     452.88: (426.75/425.73)
##     453.26: (427.13/425.73)✔
## 17: 6.5% (1054): 1 ph(s) : (1/1)--t: 17
##     31.06: (2.52/-0.16)✔
## 18: 6.8% (1066): 1 ph(s) : (1/1)--t: 18
##     29.52: (-1.54/-2.21)✔
## 19: 7.2% (1067): 1 ph(s) : (1/1)--t: 19
##     29.5: (-0.02/-0.02)✔
## 20: 7.6% (1040): 2 ph(s) : (1/2)--t: 20
##     30.86: (1.35/0.91)✔
## 21: 8% (1051): 2 ph(s) : (1/2)--t: 21
##     32.35: (1.49/-0.02)✔
## 22: 8.4% (1056): 1 ph(s) : (1/1)--t: 22
##     32.38: (0.03/-0.87)✔
## 23: 8.7% (1053): 3 ph(s) : (1/3)--t: 23
##     35.15: (2.77/2.94)✔
## 24: 9.1% (1060): 1 ph(s) : (1/1)--t: 24
##     36.58: (1.43/-1.02)✔
## 25: 9.5% (1041): 2 ph(s) : (1/2)--t: 25
##     37.64: (1.05/-0.7)✔
## 26: 9.9% (1042): 1 ph(s) : (1/1)--t: 26
##     36.69: (-0.94/-1.21)✔
## 27: 10.3% (1029): 1 ph(s) : (1/1)--t: 27
##     39.05: (2.36/1.14)✔
## 28: 10.6% (1045): 1 ph(s) : (1/1)--t: 28
##     39.5: (0.45/-0.77)✔
## 29: 11% (1063): 1 ph(s) : (1/1)--t: 29
##     39.59: (0.09/-0.25)✔
## 30: 11.4% (1057): 1 ph(s) : (1/1)--t: 30
##     40.07: (0.48/0.32)✔
## 31: 11.8% (1049): 1 ph(s) : (1/1)--t: 31
##     40.55: (0.48/0)✔
## 32: 12.2% (1048): 1 ph(s) : (1/1)--t: 32
##     40.58: (0.03/-0.45)✔
## 33: 12.5% (1047): 1 ph(s) : (1/1)--t: 33
##     40.61: (0.03/0)✔
## 34: 12.9% (1050): 1 ph(s) : (1/1)--t: 34
##     40.64: (0.03/0)✔
## 35: 13.3% (1035): 1 ph(s) : (1/1)--t: 35
##     42.2: (1.55/1.26)✔
## 36: 13.7% (1026): 1 ph(s) : (1/1)--t: 36
##     43.83: (1.64/-0.85)✔
## 37: 14.1% (1030): 1 ph(s) : (1/1)--t: 37
##     45.28: (1.45/0.7)✔
## 38: 14.4% (1038): 1 ph(s) : (1/1)--t: 38
##     46.41: (1.13/-0.2)✔
## 39: 14.8% (1032): 1 ph(s) : (1/1)--t: 39
##     46.39: (-0.02/-0.92)✔
## 40: 15.2% (1033): 1 ph(s) : (1/1)--t: 40
##     46.4: (0.01/0)✔
## 41: 15.6% (1036): 1 ph(s) : (1/1)--t: 41
##     46.68: (0.28/0.14)✔
## 42: 16% (1043): 1 ph(s) : (1/1)--t: 42
##     46.67: (-0.01/-0.11)✔
## 43: 16.3% (1037): 1 ph(s) : (1/1)--t: 43
##     45.61: (-1.06/0)✔
## 44: 16.7% (1031): 1 ph(s) : (1/1)--t: 44
##     46.31: (0.7/0.25)✔
## 45: 17.1% (1020): 1 ph(s) : (1/1)--t: 45
##     47.98: (1.67/1.51)✔
## 46: 17.5% (1028): 1 ph(s) : (1/1)--t: 46
##     47.8: (-0.18/-1.68)✔
## 47: 17.9% (1021): 1 ph(s) : (1/1)--t: 47
##     48.86: (1.06/0.35)✔
## 48: 18.3% (1039): 1 ph(s) : (1/1)--t: 48
##     48.6: (-0.25/-0.26)✔
## 49: 18.6% (1055): 3 ph(s) : (1/3)--t: 49
##     77.26: (28.66/24.49)
##     80.61: (32/26.86)
##     82.96: (34.36/29.32)✖
## 49: not included (map expansion)
## 50: 19% (1024): 1 ph(s) : (1/1)--t: 49
##     48.72: (0.12/-0.08)✔
## 51: 19.4% (1027): 1 ph(s) : (1/1)--t: 50
##     48.75: (0.03/-0.08)✔
## 52: 19.8% (1003): 6 ph(s) : (1/6)--t: 51
##     53.45: (4.7/3.4)✔
## 53: 20.2% (1015): 1 ph(s) : (1/1)--t: 52
##     53.77: (0.32/-3.11)✔
## 54: 20.5% (1012): 1 ph(s) : (1/1)--t: 53
##     53.98: (0.21/-0.26)✔
## 55: 20.9% (1014): 1 ph(s) : (1/1)--t: 54
##     55.53: (1.55/0.99)✔
## 56: 21.3% (1002): 1 ph(s) : (1/1)--t: 55
##     56.5: (0.97/0.2)✔
## 57: 21.7% (997): 3 ph(s) : (1/3)--t: 56
##     60.79: (4.3/2.59)✔
## 58: 22.1% (989): 1 ph(s) : (1/1)--t: 57
##     62.29: (1.5/-1.67)✔
## 59: 22.4% (1023): 1 ph(s) : (1/1)--t: 58
##     64.59: (2.3/-0.07)✔
## 60: 22.8% (1000): 1 ph(s) : (1/1)--t: 59
##     66.98: (2.38/0.24)✔
## 61: 23.2% (996): 1 ph(s) : (1/1)--t: 60
##     67.35: (0.38/-1.41)✔
## 62: 23.6% (1013): 4 ph(s) : (1/4)--t: 61
##     69.82: (2.47/0.77)✔
## 63: 24% (983): 1 ph(s) : (1/1)--t: 62
##     72.32: (2.5/0.47)✔
## 64: 24.3% (1018): 1 ph(s) : (1/1)--t: 63
##     72.58: (0.25/-2.01)✔
## 65: 24.7% (1005): 1 ph(s) : (1/1)--t: 64
##     73.51: (0.93/0.76)✔
## 66: 25.1% (1007): 1 ph(s) : (1/1)--t: 65
##     73.99: (0.49/-0.44)✔
## 67: 25.5% (1006): 1 ph(s) : (1/1)--t: 66
##     74.02: (0.02/-0.45)✔
## 68: 25.9% (1017): 1 ph(s) : (1/1)--t: 67
##     74.05: (0.03/0)✔
## 69: 26.2% (1009): 1 ph(s) : (1/1)--t: 68
##     74.08: (0.03/0)✔
## 70: 26.6% (992): 1 ph(s) : (1/1)--t: 69
##     75.43: (1.35/0.18)✔
## 71: 27% (1004): 1 ph(s) : (1/1)--t: 70
##     75.73: (0.3/-0.01)✔
## 72: 27.4% (988): 1 ph(s) : (1/1)--t: 71
##     76.4: (0.67/0.09)✔
## 73: 27.8% (998): 1 ph(s) : (1/1)--t: 72
##     79.45: (3.05/0.99)✔
## 74: 28.1% (990): 1 ph(s) : (1/1)--t: 73
##     79.46: (0.01/-1.22)✔
## 75: 28.5% (1001): 1 ph(s) : (1/1)--t: 74
##     79.84: (0.38/0.64)✔
## 76: 28.9% (985): 1 ph(s) : (1/1)--t: 75
##     82.32: (2.48/1.49)✔
## 77: 29.3% (994): 1 ph(s) : (1/1)--t: 76
##     81.85: (-0.47/-2.16)✔
## 78: 29.7% (993): 1 ph(s) : (1/1)--t: 77
##     81.88: (0.03/-0.01)✔
## 79: 30% (1016): 2 ph(s) : (1/2)--t: 78
##     104.71: (22.83/17.96)
##     113.64: (31.76/29.19)✖
## 79: not included (map expansion)
## 80: 30.4% (979): 1 ph(s) : (1/1)--t: 78
##     84.86: (2.98/1.65)✔
## 81: 30.8% (987): 1 ph(s) : (1/1)--t: 79
##     85.32: (0.46/-1.35)✔
## 82: 31.2% (977): 1 ph(s) : (1/1)--t: 80
##     85.68: (0.35/0.17)✔
## 83: 31.6% (974): 1 ph(s) : (1/1)--t: 81
##     85.61: (-0.07/-0.49)✔
## 84: 31.9% (980): 1 ph(s) : (1/1)--t: 82
##     85.97: (0.36/0.86)✔
## 85: 32.3% (991): 1 ph(s) : (1/1)--t: 83
##     86.81: (0.84/-0.09)✔
## 86: 32.7% (984): 1 ph(s) : (1/1)--t: 84
##     115.93: (29.12/24.95)✖
## 86: not included (map expansion)
## 87: 33.1% (976): 1 ph(s) : (1/1)--t: 84
##     89.35: (2.54/1.64)✔
## 88: 33.5% (973): 1 ph(s) : (1/1)--t: 85
##     89.61: (0.26/-2.36)✔
## 89: 33.8% (986): 2 ph(s) : (1/2)--t: 86
##     113.48: (23.87/18.01)
##     112.23: (22.62/17.1)✖
## 89: not included (map expansion)
## 90: 34.2% (969): 1 ph(s) : (1/1)--t: 86
##     90.1: (0.49/0.43)✔
## 91: 34.6% (961): 1 ph(s) : (1/1)--t: 87
##     89.93: (-0.17/-0.26)✔
## 92: 35% (970): 1 ph(s) : (1/1)--t: 88
##     90.83: (0.9/0.25)✔
## 93: 35.4% (968): 1 ph(s) : (1/1)--t: 89
##     91.33: (0.49/-0.07)✔
## 94: 35.7% (931): 1 ph(s) : (1/1)--t: 90
##     92.51: (1.18/-0.09)✔
## 95: 36.1% (978): 1 ph(s) : (1/1)--t: 91
##     92.54: (0.03/-0.3)✔
## 96: 36.5% (958): 1 ph(s) : (1/1)--t: 92
##     92.57: (0.03/0)✔
## 97: 36.9% (975): 1 ph(s) : (1/1)--t: 93
##     92.64: (0.07/0.05)✔
## 98: 37.3% (953): 1 ph(s) : (1/1)--t: 94
##     95.57: (2.93/2.46)✔
## 99: 37.6% (965): 1 ph(s) : (1/1)--t: 95
##     94.93: (-0.65/-2.47)✔
## 100: 38% (962): 2 ph(s) : (1/2)--t: 96
##     95.21: (0.28/0.32)✔
## 101: 38.4% (957): 1 ph(s) : (1/1)--t: 97
##     99.47: (4.26/2.57)✔
## 102: 38.8% (963): 1 ph(s) : (1/1)--t: 98
##     98.82: (-0.65/-2.93)✔
## 103: 39.2% (951): 1 ph(s) : (1/1)--t: 99
##     99.65: (0.83/0.29)✔
## 104: 39.5% (944): 1 ph(s) : (1/1)--t: 100
##     99.68: (0.03/-0.29)✔
## 105: 39.9% (960): 1 ph(s) : (1/1)--t: 101
##     99.82: (0.14/0.14)✔
## 106: 40.3% (934): 1 ph(s) : (1/1)--t: 102
##     100.53: (0.71/0.27)✔
## 107: 40.7% (956): 1 ph(s) : (1/1)--t: 103
##     101.98: (1.45/0.14)✔
## 108: 41.1% (955): 1 ph(s) : (1/1)--t: 104
##     103.01: (1.02/0.16)✔
## 109: 41.4% (964): 1 ph(s) : (1/1)--t: 105
##     105.52: (2.51/1.82)✔
## 110: 41.8% (959): 1 ph(s) : (1/1)--t: 106
##     105.33: (-0.19/-2.21)✔
## 111: 42.2% (952): 1 ph(s) : (1/1)--t: 107
##     107.72: (2.39/1.86)✔
## 112: 42.6% (932): 1 ph(s) : (1/1)--t: 108
##     109.28: (1.56/-1.18)✔
## 113: 43% (954): 1 ph(s) : (1/1)--t: 109
##     109.03: (-0.24/-0.82)✔
## 114: 43.3% (946): 1 ph(s) : (1/1)--t: 110
##     112.74: (3.71/2.58)✔
## 115: 43.7% (933): 1 ph(s) : (1/1)--t: 111
##     114.95: (2.2/-0.58)✔
## 116: 44.1% (947): 1 ph(s) : (1/1)--t: 112
##     114.88: (-0.06/-1.56)✔
## 117: 44.5% (950): 1 ph(s) : (1/1)--t: 113
##     115.71: (0.82/0.11)✔
## 118: 44.9% (942): 1 ph(s) : (1/1)--t: 114
##     115.99: (0.28/-0.48)✔
## 119: 45.2% (940): 1 ph(s) : (1/1)--t: 115
##     115.99: (0/-0.27)✔
## 120: 45.6% (941): 1 ph(s) : (1/1)--t: 116
##     116: (0.01/0)✔
## 121: 46% (949): 1 ph(s) : (1/1)--t: 117
##     116.01: (0.01/0)✔
## 122: 46.4% (945): 1 ph(s) : (1/1)--t: 118
##     116.02: (0.01/0)✔
## 123: 46.8% (943): 1 ph(s) : (1/1)--t: 119
##     116.03: (0.01/0)✔
## 124: 47.1% (948): 1 ph(s) : (1/1)--t: 120
##     116.04: (0.01/0)✔
## 125: 47.5% (935): 1 ph(s) : (1/1)--t: 121
##     116.45: (0.41/0.1)✔
## 126: 47.9% (936): 1 ph(s) : (1/1)--t: 122
##     116.74: (0.29/0.18)✔
## 127: 48.3% (939): 1 ph(s) : (1/1)--t: 123
##     118.01: (1.27/0.55)✔
## 128: 48.7% (937): 1 ph(s) : (1/1)--t: 124
##     117.9: (-0.1/-0.82)✔
## 129: 49% (938): 1 ph(s) : (1/1)--t: 125
##     119.04: (1.13/0.88)✔
## 130: 49.4% (925): 1 ph(s) : (1/1)--t: 126
##     119.34: (0.31/-0.75)✔
## 131: 49.8% (922): 1 ph(s) : (1/1)--t: 127
##     120.8: (1.46/0.77)✔
## 132: 50.2% (929): 1 ph(s) : (1/1)--t: 128
##     121.27: (0.47/-0.08)✔
## 133: 50.6% (921): 1 ph(s) : (1/1)--t: 129
##     121.58: (0.31/-0.66)✔
## 134: 51% (927): 1 ph(s) : (1/1)--t: 130
##     122.32: (0.74/0.11)✔
## 135: 51.3% (920): 1 ph(s) : (1/1)--t: 131
##     122.69: (0.37/-0.1)✔
## 136: 51.7% (924): 1 ph(s) : (1/1)--t: 132
##     122.71: (0.02/-0.15)✔
## 137: 52.1% (926): 1 ph(s) : (1/1)--t: 133
##     123.13: (0.42/0.08)✔
## 138: 52.5% (919): 1 ph(s) : (1/1)--t: 134
##     123.51: (0.39/0.32)✔
## 139: 52.9% (918): 1 ph(s) : (1/1)--t: 135
##     124.07: (0.56/0)✔
## 140: 53.2% (928): 1 ph(s) : (1/1)--t: 136
##     124.22: (0.15/0)✔
## 141: 53.6% (915): 1 ph(s) : (1/1)--t: 137
##     126.94: (2.72/2)✔
## 142: 54% (917): 1 ph(s) : (1/1)--t: 138
##     126.48: (-0.46/-2.24)✔
## 143: 54.4% (909): 1 ph(s) : (1/1)--t: 139
##     126.32: (-0.15/-0.09)✔
## 144: 54.8% (911): 1 ph(s) : (1/1)--t: 140
##     126.36: (0.03/-0.05)✔
## 145: 55.1% (912): 1 ph(s) : (1/1)--t: 141
##     126.39: (0.03/0)✔
## 146: 55.5% (914): 1 ph(s) : (1/1)--t: 142
##     127.09: (0.7/0.33)✔
## 147: 55.9% (916): 1 ph(s) : (1/1)--t: 143
##     127.65: (0.56/0.48)✔
## 148: 56.3% (930): 1 ph(s) : (1/1)--t: 144
##     128.64: (0.99/-0.08)✔
## 149: 56.7% (898): 1 ph(s) : (1/1)--t: 145
##     131.36: (2.72/0.44)✔
## 150: 57% (899): 1 ph(s) : (1/1)--t: 146
##     134.9: (3.54/0.86)✔
## 151: 57.4% (908): 1 ph(s) : (1/1)--t: 147
##     134.44: (-0.46/-1.89)✔
## 152: 57.8% (906): 1 ph(s) : (1/1)--t: 148
##     134.46: (0.02/0.06)✔
## 153: 58.2% (891): 1 ph(s) : (1/1)--t: 149
##     134.96: (0.5/-0.17)✔
## 154: 58.6% (890): 1 ph(s) : (1/1)--t: 150
##     136: (1.03/0.44)✔
## 155: 58.9% (896): 1 ph(s) : (1/1)--t: 151
##     136.2: (0.21/-0.18)✔
## 156: 59.3% (893): 1 ph(s) : (1/1)--t: 152
##     138.41: (2.21/1.89)✔
## 157: 59.7% (907): 1 ph(s) : (1/1)--t: 153
##     137.54: (-0.87/-2.18)✔
## 158: 60.1% (894): 1 ph(s) : (1/1)--t: 154
##     138.02: (0.48/0.44)✔
## 159: 60.5% (903): 1 ph(s) : (1/1)--t: 155
##     138.05: (0.03/-0.45)✔
## 160: 60.8% (910): 1 ph(s) : (1/1)--t: 156
##     141.59: (3.54/3.13)✔
## 161: 61.2% (879): 2 ph(s) : (1/2)--t: 157
##     142.1: (0.51/-1.39)✔
## 162: 61.6% (900): 1 ph(s) : (1/1)--t: 158
##     142.71: (0.61/-0.02)✔
## 163: 62% (881): 1 ph(s) : (1/1)--t: 159
##     142.85: (0.13/-1.36)✔
## 164: 62.4% (887): 1 ph(s) : (1/1)--t: 160
##     143.3: (0.46/-0.18)✔
## 165: 62.7% (902): 1 ph(s) : (1/1)--t: 161
##     143.1: (-0.2/-0.18)✔
## 166: 63.1% (877): 1 ph(s) : (1/1)--t: 162
##     144.68: (1.57/0.68)✔
## 167: 63.5% (884): 1 ph(s) : (1/1)--t: 163
##     145.11: (0.43/-0.42)✔
## 168: 63.9% (892): 1 ph(s) : (1/1)--t: 164
##     147.4: (2.29/1.61)✔
## 169: 64.3% (880): 1 ph(s) : (1/1)--t: 165
##     148.12: (0.72/-1.16)✔
## 170: 64.6% (885): 1 ph(s) : (1/1)--t: 166
##     148.33: (0.2/-0.7)✔
## 171: 65% (883): 1 ph(s) : (1/1)--t: 167
##     148.62: (0.29/0.13)✔
## 172: 65.4% (878): 1 ph(s) : (1/1)--t: 168
##     148.92: (0.31/-0.05)✔
## 173: 65.8% (886): 1 ph(s) : (1/1)--t: 169
##     148.93: (0.01/-0.09)✔
## 174: 66.2% (888): 1 ph(s) : (1/1)--t: 170
##     149.38: (0.44/0.13)✔
## 175: 66.5% (889): 1 ph(s) : (1/1)--t: 171
##     149.08: (-0.3/-0.08)✔
## 176: 66.9% (876): 1 ph(s) : (1/1)--t: 172
##     150.14: (1.06/0.65)✔
## 177: 67.3% (874): 1 ph(s) : (1/1)--t: 173
##     151.68: (1.55/0.76)✔
## 178: 67.7% (872): 1 ph(s) : (1/1)--t: 174
##     153.34: (1.65/0.12)✔
## 179: 68.1% (871): 1 ph(s) : (1/1)--t: 175
##     152.89: (-0.45/-1.54)✔
## 180: 68.4% (869): 1 ph(s) : (1/1)--t: 176
##     152.9: (0.01/-0.04)✔
## 181: 68.8% (870): 1 ph(s) : (1/1)--t: 177
##     153.23: (0.32/0.17)✔
## 182: 69.2% (867): 1 ph(s) : (1/1)--t: 178
##     154.19: (0.97/0.81)✔
## 183: 69.6% (875): 1 ph(s) : (1/1)--t: 179
##     153.83: (-0.36/-0.92)✔
## 184: 70% (873): 1 ph(s) : (1/1)--t: 180
##     154.76: (0.93/0.65)✔
## 185: 70.3% (865): 1 ph(s) : (1/1)--t: 181
##     154.85: (0.08/-0.33)✔
## 186: 70.7% (866): 1 ph(s) : (1/1)--t: 182
##     155.14: (0.29/-0.1)✔
## 187: 71.1% (836): 1 ph(s) : (1/1)--t: 183
##     158.7: (3.55/1.83)✔
## 188: 71.5% (863): 1 ph(s) : (1/1)--t: 184
##     158.79: (0.1/-1.81)✔
## 189: 71.9% (868): 1 ph(s) : (1/1)--t: 185
##     159.23: (0.44/0.26)✔
## 190: 72.2% (859): 1 ph(s) : (1/1)--t: 186
##     160.01: (0.78/0.15)✔
## 191: 72.6% (860): 1 ph(s) : (1/1)--t: 187
##     160.04: (0.03/-0.68)✔
## 192: 73% (854): 1 ph(s) : (1/1)--t: 188
##     161.43: (1.38/0.63)✔
## 193: 73.4% (862): 1 ph(s) : (1/1)--t: 189
##     163.23: (1.81/1.25)✔
## 194: 73.8% (853): 1 ph(s) : (1/1)--t: 190
##     164.42: (1.18/-0.5)✔
## 195: 74.1% (852): 1 ph(s) : (1/1)--t: 191
##     164.43: (0.01/-1.4)✔
## 196: 74.5% (851): 1 ph(s) : (1/1)--t: 192
##     165.05: (0.62/0.62)✔
## 197: 74.9% (822): 1 ph(s) : (1/1)--t: 193
##     166.29: (1.24/0.67)✔
## 198: 75.3% (857): 1 ph(s) : (1/1)--t: 194
##     166.31: (0.02/-0.85)✔
## 199: 75.7% (850): 1 ph(s) : (1/1)--t: 195
##     167.51: (1.2/0.72)✔
## 200: 76% (856): 1 ph(s) : (1/1)--t: 196
##     168.28: (0.76/-0.3)✔
## 201: 76.4% (855): 1 ph(s) : (1/1)--t: 197
##     168.53: (0.26/-0.61)✔
## 202: 76.8% (858): 1 ph(s) : (1/1)--t: 198
##     169.08: (0.55/0.27)✔
## 203: 77.2% (826): 1 ph(s) : (1/1)--t: 199
##     169.88: (0.8/0.14)✔
## 204: 77.6% (864): 1 ph(s) : (1/1)--t: 200
##     169.86: (-0.02/0.07)✔
## 205: 77.9% (848): 1 ph(s) : (1/1)--t: 201
##     171.62: (1.76/0.53)✔
## 206: 78.3% (849): 1 ph(s) : (1/1)--t: 201
##     170.51: (-1.11/-1.28)✔
## 207: 78.7% (861): 1 ph(s) : (1/1)--t: 201
##     168.93: (-1.58/0.17)✔
## 208: 79.1% (808): 1 ph(s) : (1/1)--t: 201
##     173.06: (4.13/3.03)✔
## 209: 79.5% (811): 1 ph(s) : (1/1)--t: 201
##     162.38: (-10.67/-2.98)✔
## 210: 79.8% (847): 1 ph(s) : (1/1)--t: 201
##     163.36: (0.97/0.21)✔
## 211: 80.2% (838): 1 ph(s) : (1/1)--t: 201
##     159.21: (-4.15/-0.43)✔
## 212: 80.6% (823): 1 ph(s) : (1/1)--t: 201
##     159.32: (0.11/0.25)✔
## 213: 81% (813): 1 ph(s) : (1/1)--t: 201
##     156.37: (-2.95/-0.07)✔
## 214: 81.4% (816): 1 ph(s) : (1/1)--t: 201
##     154.61: (-1.75/-0.17)✔
## 215: 81.7% (812): 1 ph(s) : (1/1)--t: 201
##     152.46: (-2.15/-0.01)✔
## 216: 82.1% (829): 1 ph(s) : (1/1)--t: 201
##     150.32: (-2.14/0.14)✔
## 217: 82.5% (828): 1 ph(s) : (1/1)--t: 201
##     149.79: (-0.53/-0.11)✔
## 218: 82.9% (837): 1 ph(s) : (1/1)--t: 201
##     149.84: (0.05/0.29)✔
## 219: 83.3% (840): 1 ph(s) : (1/1)--t: 201
##     149.84: (0/-0.29)✔
## 220: 83.7% (825): 1 ph(s) : (1/1)--t: 201
##     149.29: (-0.55/0)✔
## 221: 84% (815): 1 ph(s) : (1/1)--t: 201
##     147.95: (-1.34/0.57)✔
## 222: 84.4% (814): 1 ph(s) : (1/1)--t: 201
##     147: (-0.95/-0.57)✔
## 223: 84.8% (846): 1 ph(s) : (1/1)--t: 201
##     147.05: (0.05/0.04)✔
## 224: 85.2% (841): 1 ph(s) : (1/1)--t: 201
##     146.69: (-0.36/-0.04)✔
## 225: 85.6% (804): 1 ph(s) : (1/1)--t: 201
##     146.19: (-0.49/0.1)✔
## 226: 85.9% (820): 3 ph(s) : (1/3)--t: 201
##     146.59: (0.39/0.18)✔
## 227: 86.3% (819): 1 ph(s) : (1/1)--t: 201
##     145.28: (-1.31/-0.06)✔
## 228: 86.7% (818): 1 ph(s) : (1/1)--t: 201
##     144.51: (-0.76/0.54)✔
## 229: 87.1% (809): 1 ph(s) : (1/1)--t: 201
##     142.73: (-1.79/-0.49)✔
## 230: 87.5% (805): 1 ph(s) : (1/1)--t: 201
##     141.29: (-1.43/-0.3)✔
## 231: 87.8% (827): 1 ph(s) : (1/1)--t: 201
##     142.88: (1.58/0.83)✔
## 232: 88.2% (803): 1 ph(s) : (1/1)--t: 201
##     142.69: (-0.19/-0.43)✔
## 233: 88.6% (824): 1 ph(s) : (1/1)--t: 201
##     142.79: (0.11/-0.21)✔
## 234: 89% (797): 1 ph(s) : (1/1)--t: 201
##     145.8: (3/1.92)✔
## 235: 89.4% (807): 1 ph(s) : (1/1)--t: 201
##     148.14: (2.34/-1.92)✔
## 236: 89.7% (806): 1 ph(s) : (1/1)--t: 201
##     148.14: (0/-0.17)✔
## 237: 90.1% (800): 1 ph(s) : (1/1)--t: 201
##     149.61: (1.47/0.79)✔
## 238: 90.5% (810): 1 ph(s) : (1/1)--t: 201
##     150.39: (0.78/-0.17)✔
## 239: 90.9% (801): 1 ph(s) : (1/1)--t: 201
##     152.41: (2.02/1.76)✔
## 240: 91.3% (796): 3 ph(s) : (1/3)--t: 201
##     151.92: (-0.49/-1.93)✔
## 241: 91.6% (788): 4 ph(s) : (1/4)--t: 201
##     154.73: (2.8/2.28)✔
## 242: 92% (789): 1 ph(s) : (1/1)--t: 201
##     153.94: (-0.78/-2.76)✔
## 243: 92.4% (787): 1 ph(s) : (1/1)--t: 201
##     153.9: (-0.04/0)✔
## 244: 92.8% (793): 1 ph(s) : (1/1)--t: 201
##     153.67: (-0.24/0)✔
## 245: 93.2% (792): 1 ph(s) : (1/1)--t: 201
##     153.67: (0.01/0)✔
## 246: 93.5% (791): 1 ph(s) : (1/1)--t: 201
##     153.23: (-0.45/0)✔
## 247: 93.9% (790): 1 ph(s) : (1/1)--t: 201
##     153.32: (0.1/0.09)✔
## 248: 94.3% (795): 1 ph(s) : (1/1)--t: 201
##     152.49: (-0.83/0.28)✔
## 249: 94.7% (794): 1 ph(s) : (1/1)--t: 201
##     150.01: (-2.48/0.73)✔
## 250: 95.1% (783): 1 ph(s) : (1/1)--t: 201
##     151.67: (1.66/1.39)✔
## 251: 95.4% (802): 1 ph(s) : (1/1)--t: 201
##     152.44: (0.76/-2.12)✔
## 252: 95.8% (774): 1 ph(s) : (1/1)--t: 201
##     154.07: (1.63/0.44)✔
## 253: 96.2% (780): 1 ph(s) : (1/1)--t: 201
##     152.35: (-1.72/-0.78)✔
## 254: 96.6% (772): 4 ph(s) : (1/4)--t: 201
##     155.13: (2.79/0.53)✔
## 255: 97% (770): 1 ph(s) : (1/1)--t: 201
##     152.6: (-2.54/-0.53)✔
## 256: 97.3% (776): 1 ph(s) : (1/1)--t: 201
##     152.7: (0.1/0.6)✔
## 257: 97.7% (779): 1 ph(s) : (1/1)--t: 201
##     152.43: (-0.26/-0.63)✔
## 258: 98.1% (765): 1 ph(s) : (1/1)--t: 201
##     152.25: (-0.18/0.4)✔
## 259: 98.5% (778): 1 ph(s) : (1/1)--t: 201
##     152.78: (0.53/0.51)✔
## 260: 98.9% (777): 1 ph(s) : (1/1)--t: 201
##     151.8: (-0.99/-0.47)✔
## 261: 99.2% (782): 1 ph(s) : (1/1)--t: 201
##     151.19: (-0.6/0.01)✔
## 262: 99.6% (786): 1 ph(s) : (1/1)--t: 201
##     151.67: (0.48/0.88)✔
## 263: 100% (769): 1 ph(s) : (1/1)--t: 201
##     149.35: (-2.33/-0.18)✔
```

```
## ═════════════════════════════════════════════════════════════ Reestimating final recombination fractions ══
```

```
## Markers in the initial sequence: 263
## Maped markers                  : 259 (98.5%)
```

```
## ═══════════════════════════════════════════════════════════════════════════════════════════════════════════
```


```r
plot(lg3.map.mds)
```

![](README_files/figure-html/lg3_map_mds_plot-1.png)<!-- -->

```r
plot_map_list(list(lg3.map, lg3.map.mds), col = c("red", "blue"), title = "")
```

![](README_files/figure-html/lg3_map_mds_plot-2.png)<!-- -->

Comparing the genomic-based and the MDS-based maps, the first one is smaller, with the same number of markers, indicating a better result. However, to formally compare both maps, we need to select the markers present in both maps should be the same. Interestingly enough, both procedures included the same markers in the final map. However, we provide the code to perform the comparison even if the maps share a subset of markers


```r
mrks.in.gen<-intersect(lg3.map$maps[[1]]$seq.num, lg3.map.mds$maps[[1]]$seq.num)
mrks.in.mds<-intersect(lg3.map.mds$maps[[1]]$seq.num, lg3.map$maps[[1]]$seq.num)
if(cor(mrks.in.gen, mrks.in.mds) < 0){
  mrks.in.mds <- rev(mrks.in.mds)
  lg3.map.mds <- rev_map(lg3.map.mds)
}
map.comp.3.gen<-get_submap(input.map = lg3.map, match(mrks.in.gen, lg3.map$maps[[1]]$seq.num), verbose = FALSE)
map.comp.3.mds<-get_submap(input.map = lg3.map.mds, match(mrks.in.mds, lg3.map.mds$maps[[1]]$seq.num), verbose = FALSE)
dist.3.gen<-extract_map(lg3.map)
dist.3.mds<-extract_map(lg3.map.mds)
names(dist.3.gen)<-map.comp.3.gen$maps[[1]]$seq.num
names(dist.3.mds)<-map.comp.3.mds$maps[[1]]$seq.num       
matplot(t(data.frame(dist.3.gen,dist.3.mds[names(dist.3.gen)])), 
        type="b", pch="_", col=1, lty=1, lwd = .5, xlab= "", 
        ylab="Marker position (cM)", axes = F)
axis(2)
mtext(text = round(map.comp.3.gen$maps[[1]]$loglike,1), side = 1, adj = 0)
mtext(text = round(map.comp.3.mds$maps[[1]]$loglike,1), side = 1, adj = 1)
mtext(text = "Genomic", side = 3, adj = 0)
mtext(text = "MDS", side = 3, adj = 1)
```

![](README_files/figure-html/map_comp-1.png)<!-- -->

Notice that we can observe the same local inversions shown in the dot plots presented earlier. In this case, the log-likelihood of the genomic order is higher than the one obtained using the MDS order. Thus for this linkage group, we opted to use the genomic-based map.

# Parallel map construction

Now we will apply the mapping procedure to all linkage groups using parallelization. We will use the genomic order as an example. However, the user should compare the MDS and genomic maps in a similar way we did in the previous example. Also, some species may not have a useful reference genome. Thus, the MDS order is a valuable alternative.  


```r
## Performing parallel computation
my.phase.func<-function(X){
  x<-est_rf_hmm_sequential(input.seq = X$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 200,
                                info.tail = TRUE, 
                                twopt = X$tpt,
                                sub.map.size.diff.limit = 8, 
                                phase.number.limit = 10,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
  return(x)
}
system.time({
  cl <- parallel::makeCluster(12)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.dose.filt")
  MAPs <- parallel::parLapply(cl,LGS,my.phase.func)
  parallel::stopCluster(cl)
})
```

We can plot a traditional linkage map including all linkage groups using the function `plot_map_list`


```r
plot_map_list(MAPs, col = "ggstyle")
```

![](README_files/figure-html/print_maps-1.png)<!-- -->

Similarly we did for linkage group 3, let us consider a global genotyping error of 5% to reestimate the final maps.


```r
my.error.func<-function(X){
  x<-est_full_hmm_with_global_error(input.map = X, 
                                    error = 0.05, 
                                    tol = 10e-4, 
                                    verbose = FALSE)
  return(x)
}
system.time({
  cl <- parallel::makeCluster(12)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.dose.filt")
  MAP.err <- parallel::parLapply(cl,MAPs,my.error.func)
  parallel::stopCluster(cl)
})
```

```
##    user  system elapsed 
##   0.182   2.107 147.633
```

Let us compare both results


```r
all.MAPs <- NULL
for(i in 1:12) 
  all.MAPs<-c(all.MAPs, MAPs[i], MAP.err[i])
plot_map_list(map.list = all.MAPs, col = rep(c("#E69F00", "#56B4E9"), 12))
```

![](README_files/figure-html/final_plot-1.png)<!-- -->

As our final map, we will chose the the one that includes the modeling of genotype errors.


```r
plot_map_list(MAP.err)
```

![](README_files/figure-html/print_err_maps-1.png)<!-- -->

# Genotype conditional probabilities

In order to use the genetic map in [QTLpoly](https://github.com/guilherme-pereira/QTLpoly), we need to obtain the conditional probability of all possible 36 genotypes along the 12 linkage groups for all individuals in the full-sib population. Again, let us use the function `calc_genoprob_error`. 


```r
genoprob.err <- vector("list", 12)
for(i in 1:12)
genoprob.err[[i]] <- calc_genoprob_error(input.map = MAP.err[[i]], error = 0.05)
```

Let us display the results for all linkage groups, individual 1


```r
ind <- 1
dg <- sapply(genoprob.err, function (x) max(x$map))
dg <- dg/max(dg)
layout(matrix(1:12, ncol = 4), widths = dg)
for(i in 1:12)
{
  d <- genoprob.err[[i]]$map
  image(t(genoprob.err[[i]]$probs[,,ind]),
        col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
        axes=FALSE,
        xlab = "Markers",
        ylab = "",
        main = paste("LG", i))
  axis(side = 1, at = d/max(d),
       labels =rep("", length(d)), las=2)
}
```

![](README_files/figure-html/plot_genoprob_all-1.png)<!-- -->

With the conditional probabilities computed, it is possible to use the object `genoprob.err` alongside with phenotypic data as the input of the software [QTLpoly](https://github.com/guilherme-pereira/QTLpoly).



