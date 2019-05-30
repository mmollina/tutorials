---
title: "Building a genetic map of an autotetraploid population using MAPpoly"
author: "Marcelo Mollinari, Guilhereme Pereira, A Augusto F Garcia and Zhao-Bang Zeng"
date: '2019-05-06'
output:
  html_document:
    highlight: pygments
    keep_md: yes
    theme: flatly
    toc: yes
    toc_depth: '2'
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

Here we assume that the genotypic data is available and in the format required by `mappoly`. The primary purpose of this tutorial is to show some functions available in `mappoly` and how to use them sequentially to construct a genetic map. The derivation of the HMM used in `mappoly` is presented in [@Mollinari2018](https://doi.org/10.1101/415232 ).

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

# Loading, inspecting and filtering data

In this tutorial, we will construct a genetic map of the B2721 population which is a cross between two tetraploid potato varieties: Atlantic $\times$ B1829-5. The population comprises 160 offsprings genotyped with the SolCAP Infinium 8303 potato array. The data set also contains the genomic order of the SNPs from the _Solanum tuberosum_ genome version 4.03. The genotype calling was performed using fitTetra R package using [this pipeline](https://github.com/mmollina/Autopolyploid_Linkage/blob/master/src/solcap_map_construction/snp_calling/genotype_calling_public_data_fittetra.R). Another option would be to use ClusterCall and [this pipeline](https://mmollina.github.io/tutorials/solcap/solcap_example.html). MAPpoly accepts two types of input files: the first one can be imported using the function `read_geno` contains and a header and a matrix of genotypes for each marker and offspring. A detailed description of the input file can be found in the function documentation (`?read_geno`). The second type of file contains the same header and a data frame with the (posterior) probability distribution for each combination of marker and offspring. It can be imported using the function `read_geno_dist`. Figure 1 shows an example of a data set containing genotype probabilities.

<div class="figure" style="text-align: center">
<img src="image/dataset_tetra.png" alt="Figure 1: Example of data set containing genotype probabilities" width="40%" />
<p class="caption">Figure 1: Example of data set containing genotype probabilities</p>
</div>

## Loading the data set

Now, let us read the data set using the function `read_geno_dist`.


```r
#(~46 seconds)
dat.pot<-read_geno_dist(file.in = "SolCAP", prob.thres = 0.9)
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

The loaded data set contains 4017 markers and 160 individuals. Notice that the program performed some filtering, eliminating non-informative markers (monomorphic) and also computed the p-value associated to the chi-square test for the expected segregation patterns under Mendelian inheritance for all markers. To inspect the data set, we just type


```r
print(dat.pot, detailed = TRUE)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             4017 
##     Missing data under 0.9 prob. threshold: 2.17%
## 
##     ----------
##     No. markers per sequence:
##        seq No.mrk
##          1    411
##          2    281
##          3    318
##          4    432
##          5    291
##          6    433
##          7    409
##          8    279
##          9    310
##         10    245
##         11    284
##         12    221
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
The program prints a summary of the data set showing the ploidy level, the number of individuals, the number of markers and the percentage of missing data under a probability threshold (2.17%). In this case, data points with a maximum probability smaller than 0.9 were classified as "missing". It is important to mention that for hexaploids, you should use a lower threshold to include more multiple-dose markers in the analysis. In the next section, the number of markers per sequence (in our case, chromosomes) is shown followed by the number of markers per each dosage combination in both parents. Now, let us inspect the data set using some graphics


```r
plot(dat.pot)
```

![](README_files/figure-html/plot_datset-1.png)<!-- -->

In the bar plot in the left-hand side, is a graphical representation of the number of markers per each dosage combination in both parents. The numbers separated by a dash indicate the dose in parents $P1$ and $P2$ respectively.  The $\log_{10}(p.value)$ from a chi-square test for the expected segregation patterns under Mendelian inheritance for all markers are shown in the upper-right panel. The lower-right panel shows a graphical representation of the dosage distribution in the full-sib population. 

## Filtering

After reading the data set, let us filter out markers with more than a certain percentage of missing data using the function `filter_missing`. In the following case, we are filtering out markers with more than 20% of missing data. You can make this number bigger, say `filter.thres = 0.3` or even `0.4`, but be aware that during the linkage analysis you will need to be extra cautious about false positives and marker ordering.


```r
dat.pot.filt<-filter_missing(dat.pot, filter.thres = 0.2, inter = FALSE)
print(dat.pot.filt)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:                            4 
##     No. individuals:                         160 
##     No. markers:                             4001 
##     Missing data under 0.9 prob. threshold: 2.09%
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage combination in both parents:
##     P1 P2 freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  432
##      1  1  308
##      1  2  256
##      1  3  125
##      1  4   12
##      2  0  171
##      2  1  209
##      2  2  336
##      2  3  196
##      2  4   50
##      3  0   41
##      3  1  108
##      3  2  201
##      3  3  200
##      3  4  333
##      4  1    9
##      4  2   86
##      4  3  307
```

```r
plot(dat.pot.filt) #based in a bonferroni aproximation
```

![](README_files/figure-html/filtering_missing-1.png)<!-- -->

Notice that only 16 markers were eliminated, indicating excellent data quality. Now, let us select markers with associated p-value smaller then $0.05/4001 = 1.25 \times 10^{-5}$ (approximated Bonferroni correction)


```r
seq.filt<-filter_segregation(dat.pot.filt, chisq.pval.thres = 0.05/dat.pot.filt$n.mrk, inter = FALSE)
seq.filt<-make_seq_mappoly(seq.filt)
print(seq.filt)
```

```
## This is an object of class 'mappoly.sequence'
##     ------------------------
##     Parameters not estimated
##     ------------------------
##     Ploidy level:       4 
##     No. individuals:    160 
##     No. markers:        4001 
## 
##     ----------
##     No. markers per sequence:
##   sequence No.mrk
##          1    411
##          2    280
##          3    317
##          4    430
##          5    289
##          6    433
##          7    408
##          8    278
##          9    309
##         10    244
##         11    282
##         12    219
## 
##     ----------
##     No. of markers per dosage in both parents:
##     dP dQ freq
##      0  1  427
##      0  2  170
##      0  3   24
##      1  0  432
##      1  1  308
##      1  2  256
##      1  3  125
##      1  4   12
##      2  0  171
##      2  1  209
##      2  2  336
##      2  3  196
##      2  4   50
##      3  0   41
##      3  1  108
##      3  2  201
##      3  3  200
##      3  4  333
##      4  1    9
##      4  2   86
##      4  3  307
```

In this case, the filtering did not affect since none of the markers were below the threshold. 

## Inspecting a specific marker

Now, let us assess the information of a specific marker, in this case marker `"solcap_snp_c1_13686"`, using the function `plot_mrk_info`


```r
plot_mrk_info(input.data = dat.pot.filt, mrk = "solcap_snp_c1_13686")
```

<img src="README_files/figure-html/plot_mrk_info-1.png" width="1000px" />
In this Figure, we have the position of the markers in the data set (238), de dosage in both parents (2 and 3) the amount of missing data (1,2%), The p-value for the segregation test, the chromosome and position the SNP is located in the reference genome and finally the threshold we used to assume missing data. 

## Eliminating redundant markers

Frequently closely linked markers carry the same information about a specific genomic region and can be excluded from the analysis without modifying the results. To identify those markers, we use the function `elim_redundant`.


```r
## Filtering out redundant markers
seq.uni<-elim_redundant(input.seq = seq.filt)
plot(seq.uni)
```

![](README_files/figure-html/redundant_elimi-1.png)<!-- -->

In this case, 347 markers (9%) were eliminated. Now we select the 91% non-redundant markers using the function `make_seq_mappoly`. 


```r
## Filtering out redundant markers
seq.all.lgs<-make_seq_mappoly(seq.uni)
```
OBS: To select all markers in the data set, we use the function `make_seq_mappoly` with `arg = 'all'`. It is also possible to load data only for a specific sequence using `arg = 'seqx'`, where `x` is the number of the sequence, provided in the input file. This feature can be useful in cases where the information about chromosomes or scaffold is available. It is also possible to load specific markers using a vector of numbers indicating the positions of the markers in the data set. 

# Two-point analysis

Once the markers where selected, we need to compute the pairwise recombination fraction between all these markers (two-point analysis). First, let us load the genotype counts ($\zeta_{\mbox{T}_{k},\mbox{T}_{k^{\prime}}}(l_{P}, l_{Q})$) defined in equation 20 in [@Mollinari2018](https://doi.org/10.1101/415232 ). This object is fundamental to perform the dimension reduction of the transition space.


```r
counts<-cache_counts_twopt(input.seq = seq.all.lgs, get.from.web = TRUE)
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
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.all.lgs, 
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
## 2-2   0.00000 0.02196982 5.704851e+01
## 2-0 -45.84217 0.07569726 1.120634e+01
## 0-2 -45.84217 0.07569726 1.120634e+01
## 2-1 -56.59683 0.20434970 9.939297e+00
## 1-2 -56.59683 0.20434970 9.939297e+00
## 0-0 -57.05404 0.49995416 5.534255e-03
## 1-0 -66.53882 0.49995416 2.687607e-03
## 0-1 -66.53882 0.49995416 2.687607e-03
## 1-1 -70.86902 0.49995416 9.215219e-08
```

```r
plot(all.rf.pairwise, first.mrk = 93, second.mrk = 98)
```

<img src="README_files/figure-html/plot_twopt_example-1.png" width="500px" style="display: block; margin: auto;" />

In this case, `93-98` represents the position of the markers in the filtered data set. The name of the rows have the form `x-y`, where `x` and `y` indicate how many homologous chromosomes share the same allelic variant in parents $P1$ and $P2$, respectively (see [@Mollinari2018](https://doi.org/10.1101/415232 ) for notation). The first column indicates the LOD Score in relation to the most likely linkage phase configuration. The second column shows the recombination fraction, and the third indicates the LOD Score comparing the likelihood under no linkage ($r = 0.5$) and the estimated recombination fraction (evidence of linkage).

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
id<-get_genomic_order(seq.all.lgs)
plot(mat, ord = rownames(id), index = FALSE)
```

As expected, we can observe the block-diagonal and monotonic patterns. In the previous case, the thresholds allowed to plot almost all points in the recombination fraction matrix. The empty cells in the matrix indicate markers where it is impossible to detect recombinant events using two-point estimates (e.g., between $1 \times 0$ and $0 \times 1$ marker). Yet, if the thresholds become more stringent (higher LODs and lower rf), the matrix becomes more sparse.  

# Assembling linkage groups

The function `group_mappoly` assign markers to linkage groups using the recombination fraction matrix obtained above. The user can provide an expected number of groups or run the interactive version of the function using `inter = TRUE`. Since in this data set we expect 12 linkage groups (basic chromosome number in potato), we use `expected.groups = 12`. If the data set provides the chromosome where the markers are located, the function allows comparing the groups obtained using the pairwise recombination fraction and the chromosome information provided using the `comp.mat = TRUE`.


```r
grs <- group_mappoly(input.mat = mat,
                     input.seq = seq.all.lgs,
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
##     - Number of markers:          3654 
##     - Number of linkage groups:   12 
##     - Number of markers per linkage groups: 
##     group n.mrk
##         1   427
##         2   289
##         3   284
##         4   385
##         5   278
##         6   390
##         7   356
##         8   271
##         9   291
##        10   221
##        11   248
##        12   214
##   ------------------------------------------
##      1   2   3   4   5   6   7   8   9  10  11  12 NH
## 1  237   0   0   1   0   0   2   0   0   1   4   0  3
## 2    2 271   1   1   0   1   0   0   1   1   1   2 10
## 3    0   0 329   2   5   0   2   0   6   1   0   0 11
## 4    0   0   1 257   1   1   2   0   2   2   0   0 12
## 5    2   1   2   1 363   2   1   0   1   1   0   1 10
## 6    2   1   4   5   5 380   0   9   2   3   0   3 13
## 7    0   1   1   0   0   1 254   0   1   4   0   0  9
## 8    3   0   2   3   6   1   0 256   0   4   4   6  4
## 9    0   1   0   1   0   3   0   0 207   0   3   0  6
## 10   1   3   2   2   2   0   0   0   2 369   0   0  9
## 11   2   2   1   0   0   0   1   1   0   4 201   0  2
## 12   1   0   0   0   0   1   0   1   1   1   1 277  1
##   ------------------------------------------
```

Here, we have the 3654 markers distributed in 12 linkage groups. The rows indicate linkage groups obtained using linkage information and the columns are the chromosomes in the reference genome. Notice the diagonal indicating the concordance between the two sources of information. Now, we can plot the resulting marker cluster analysis.


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
  nm<-names(which(temp2[,1]==names(which.max(table(temp2[,1])))))
  lgtemp<-make_seq_mappoly(dat.pot.filt, nm)
  LGS[[j]]<-list(lg = lgtemp, 
                 tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}
```

Now, let us print the recombination fraction matrices or each linkage group.

![](README_files/figure-html/all_mat_rf-1.png)<!-- -->

# Estimating the map for a given order

In this section, we will use the marker order provided by the _Solanum tuberosum_ genome version 4.03. The estimation of the genetic map for a given order involves the computation of recombination fraction between adjacent markers and also finding the linkage phase configuration of those markers in the parents. The core function to perform these tasks in `mappoly` is `est_rf_hmm_sequential`. This function uses the pairwise recombination fraction as the first source of information to sequentially position allelic variants in specific homologs. For situations where pairwise analysis has limited power, the algorithm relies on the likelihood obtained through a hidden Markov model (HMM) [@Mollinari2018].  Once all markers are positioned, the final map is reconstructed using the HMM multipoint algorithm. 

![Example of linkage phase configuration estimation using sequential search space
reduction and HMM evaluation.](phasing.png)

Several arguments are available to control the inclusion and phasing of the markers in the chain.  `thres.twopt` receives the threshold to whether when the linkage phases compared via two-point analysis should be considered, and the HMM analysis should not be used to infer the linkage phase (A. K. A. $\eta$ in [@Mollinari2018](https://doi.org/10.1101/415232 )). `thres.hmm` receives the threshold for keeping competing maps computed using HMM (if the two-point analysis was not enough) in the next round of marker insertion. `extend.tail` indicates the number of markers that should be considered at the end of the chain to insert a new marker. `tol` and `tol.final` receive the desired accuracy to estimate the sub-maps during the sequential phasing procedure and the desired accuracy in the final map. `phase.number.limit` receives the limit number of linkage phase configurations to be tested using HMM. `info.tail` is a logical argument and if `TRUE` it uses the complete informative tail (last markers in the chain that allow all homologous to be distinguished in the parents) of the chain to calculate the likelihood of the linkage phases. 

First, as an example, let us estimate the map for linkage group 12. The values used in the function arguments, were obtained using a balance of processing speed and accuracy of the algorithm. As an exercise, it is interesting to try different values and check out the results. For now, let us stick with these values.


```r
# (~74 seconds)
lg12.map<-est_rf_hmm_sequential(input.seq = LGS[[12]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 20,
                                info.tail = TRUE, 
                                twopt = LGS[[12]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```

Now, we can display the results using the function `print` and `plot`. 


```r
print(lg12.map)
```

```
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 197 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 172.5
##        log-likelihood:	 -5516.29
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
```

```r
plot(lg12.map)
```

![](README_files/figure-html/map_lg12_plot-1.png)<!-- -->

## Reestimaing the map considering the genotype probability distribution 

Since our data set contains the probability distribution for each combination of marker and offspring, we will use this information to reestimate the genetic map of linkage group 12. The use of this prior information will update the genotype of the markers based on a global chromosome structure.


```r
lg12.map.prior<-reest_map_geno_dist(lg12.map, dat.dist = dat.pot.filt, verbose = FALSE, tol = 10e-4)
plot(lg12.map.prior$reestimated.map)
```

![](README_files/figure-html/map_err_lg12-1.png)<!-- -->
Notice that the resulting map is smaller than the previous and some markers were "attracted" and some markers were "repealed". We will see the effect of these changes in a future section.

## Using a global error

This is similar to the strategy of using probability distribution of the markers. In this case, instead of using a specific probability distribution as the prior information, we will use a global error of 1%


```r
lg12.map.error<-est_full_hmm_with_global_error(input.map = lg12.map, error = 0.01)
plot(lg12.map.error)
```

![](README_files/figure-html/map_prior_lg12-1.png)<!-- -->

## Posterior genotype probabilities

Now, let us compute the posterior genotype probabilities. First, let us use the function  `calc_genoprob`. This function computes the genotype probabilities using the dosage information with no probability distribution associated. Actually, in this case, $\Pr(\text{assigned probabity}) = 1.0$ and zero for the remaining classes. Thus, the genotype probability on the called marker will not be updated, since we are using a full informative prior for each marker.

using functions `calc_genoprob` and `calc_genoprob_dist`. 


```r
genoprob.lg12<-calc_genoprob(input.map = lg12.map)
```

```
## 	Ploidy level: 4
## 	Number of markers: 197
## 	Number of individuals: 160
## 	..................................................
## 	..................................................
## 	..................................................
## 	..........
```

Let us display the results for offspring 1


```r
ind <- 1
d <- genoprob.lg12$map
pr<-genoprob.lg12$probs[,,ind]
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

<img src="README_files/figure-html/plot_genoprob-1.png" width="200%" style="display: block; margin: auto;" />



Now, let us use the function `calc_genoprob_dist`, whcih allows the inclusion of genotype probability distribution 


```r
genoprob.lg12.prior<-calc_genoprob_dist(input.map = lg12.map, dat.dist = dat.pot.filt)
```

```
## Ploidy level:4
## Number of individuals:160
## 	..................................................
## 	..................................................
## 	..................................................
## 	..........
```

Again, let us display the results for offspring 1


```r
ind <- 1
d <- genoprob.lg12.prior$map
pr<-genoprob.lg12.prior$probs[,,ind]
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


Notice that ...


# Ordering markers

So far we reestimated the map using the genomic order. In real situations, unless a genomic information is provided, we need to order the markers using an optimization technique. Here, we use the MDS (multidimensional scaling) algorithm, proposed in the context of genetic mapping by [@Preedy2016]. It requires a recombination fraction matrix, which will be transformed in distance using a mapping function (in this case we use Haldane's mapping function). First, let us gather the pairwise recombination fractions for all three linkage groups


```r
mt <- lapply(LGS, function(x) rf_list_to_matrix(x$tpt))
```

Now, for each matrix contained in the object in `mt`, we use the MDS algorithm


```r
mds.ord <- lapply(mt, mds_mappoly)
```

Now, let us compare the estimated and the simulated order using 

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

Although we can observe several local inconsistencies, the global diagonal patterns indicate a very good order for all linkage groups.

Now, let us build the genetic map of linkage group 12 using the MDS order


```r
# (~74 seconds)
lg12.map.mds<-est_rf_hmm_sequential(input.seq = LGS.mds[[12]]$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 20,
                                info.tail = TRUE, 
                                twopt = LGS.mds[[12]]$tpt,
                                sub.map.size.diff.limit = 10, 
                                phase.number.limit = 20,
                                reestimate.single.ph.configuration = TRUE,
                                tol = 10e-3,
                                tol.final = 10e-4)
```


```r
plot(lg12.map.mds)
```

![](README_files/figure-html/lg12_map_mds_plot-1.png)<!-- -->

```r
print(lg12.map)
print(lg12.map.mds)
```


```r
#(intersection of 190 markers, out of 197 from genomeic order and 192 f )
mrks.in.gen<-intersect(lg12.map$maps[[1]]$seq.num, lg12.map.mds$maps[[1]]$seq.num)
mrks.in.mds<-intersect(lg12.map.mds$maps[[1]]$seq.num, lg12.map$maps[[1]]$seq.num)
map.comp.12.gen<-get_submap(input.map = lg12.map, match(mrks.in.gen, lg12.map$maps[[1]]$seq.num))
map.comp.12.mds<-get_submap(input.map = lg12.map.mds, match(mrks.in.mds, lg12.map.mds$maps[[1]]$seq.num))
dist.12.gen<-cumsum(c(0, imf_h(map.comp.12.gen$maps[[1]]$seq.rf)))
dist.12.mds<-cumsum(c(0, imf_h(map.comp.12.mds$maps[[1]]$seq.rf)))
names(dist.12.gen)<-map.comp.12.gen$maps[[1]]$seq.num
names(dist.12.mds)<-map.comp.12.mds$maps[[1]]$seq.num       
matplot(t(data.frame(dist.12.gen,dist.12.mds[names(dist.12.gen)])), 
        type="b", pch="_", col=1, lty=1, lwd = .5, xlab= "", 
        ylab="Marker position (cM)", axes = F)
axis(2)
mtext(text = round(map.comp.12.gen$maps[[1]]$loglike,1), side = 1, adj = 0)
mtext(text = round(map.comp.12.mds$maps[[1]]$loglike,1), side = 1, adj = 1)
mtext(text = "Genomic", side = 3, adj = 0)
mtext(text = "MDS", side = 3, adj = 1)
```

![](README_files/figure-html/map_comp-1.png)<!-- -->


# Parallel map construction


In the fallowing example, we use the package `paralell` to perform the construction of the three linkage groups simultaneously.


```r
## Performing parallel computation
my.phase.func<-function(X){
  x<-est_rf_hmm_sequential(input.seq = X$lg,
                                start.set = 10,
                                thres.twopt = 10, 
                                thres.hmm = 10,
                                extend.tail = 50,
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
  parallel::clusterExport(cl, "dat.pot.filt")
  MAPs <- parallel::parLapply(cl,LGS,my.phase.func)
  parallel::stopCluster(cl)
})
```



```r
MAPs.err <- MAPs.prior <- vector("list", 12)
for(i in 1:12){
  cat("LG:", i ,"\n")
   MAPs.err[[i]] <- est_full_hmm_with_global_error(input.map = MAPs[[i]], error = 0.05, tol = 10e-4, verbose = FALSE)
   MAPs.prior[[i]] <- reest_map_geno_dist(input.map = MAPs[[i]], dat.dist = dat.pot.filt, verbose = FALSE, tol = 10e-4)
}
```



```r
 sapply(MAPs.prior, function(x) plot(x$reestimated.map))
```

![](README_files/figure-html/print_maps-1.png)<!-- -->![](README_files/figure-html/print_maps-2.png)<!-- -->![](README_files/figure-html/print_maps-3.png)<!-- -->![](README_files/figure-html/print_maps-4.png)<!-- -->![](README_files/figure-html/print_maps-5.png)<!-- -->![](README_files/figure-html/print_maps-6.png)<!-- -->![](README_files/figure-html/print_maps-7.png)<!-- -->![](README_files/figure-html/print_maps-8.png)<!-- -->![](README_files/figure-html/print_maps-9.png)<!-- -->![](README_files/figure-html/print_maps-10.png)<!-- -->![](README_files/figure-html/print_maps-11.png)<!-- -->![](README_files/figure-html/print_maps-12.png)<!-- -->

The results were stored in a list format in the object `maps.given.order`. A graphical representation of the linkage groups including the linkage phase configurations and the distance between markers can be obtained using

In these figures, the red and white rectangles indicate the two possible allelic variants. Each horizontal line
containing these rectangles indicates one of the six homologous chromosomes which are grouped in homology groups. 
The top homology group represents parent $P$ and the bottom represents parent $Q$. 


###Reestimating the genetic map

Now, given the estimated order, we reestimate the final map using the function `est_rf_hmm_sequential`


```r
## Performing parallel computation
system.time({
  cl <- parallel::makeCluster(12)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "dat.pot.filt")
  MAPS.mds <- parallel::parLapply(cl,LGS.mds,my.phase.func)
  parallel::stopCluster(cl)
})
```



```r
 MAPS.mds 
```

```
## [[1]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 368 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 367.68
##        log-likelihood:	 -9865.32
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[2]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 236 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 296.02
##        log-likelihood:	 -7254.02
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[3]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 271 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 271.1
##        log-likelihood:	 -7753.35
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[4]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 348 
##     No. linkage phases:	 2 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  2
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 312.86
##        log-likelihood:	 -9089.15
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
##     Linkage phase configuration:  2
##        map length:	 319.71
##        log-likelihood:	 -9118.68
##        LOD:		 -30
##     ~~~~~~~~~~~~~~~~~~
## 
## [[5]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 244 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 207.73
##        log-likelihood:	 -6442.27
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[6]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 355 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 288.87
##        log-likelihood:	 -8325.02
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[7]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 312 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 236.7
##        log-likelihood:	 -7091.11
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[8]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 243 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 193.7
##        log-likelihood:	 -5892.63
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[9]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 256 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 239.1
##        log-likelihood:	 -6968.9
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[10]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 188 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 208.45
##        log-likelihood:	 -5786.13
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[11]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 229 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 202.46
##        log-likelihood:	 -6257.13
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[12]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 185 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 202.16
##        log-likelihood:	 -5719.35
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
```
Graphical representations


```r
sapply(MAPS.mds, plot)
```

![](README_files/figure-html/plot_mds_hmm_map-1.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-2.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-3.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-4.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-5.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-6.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-7.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-8.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-9.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-10.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-11.png)<!-- -->![](README_files/figure-html/plot_mds_hmm_map-12.png)<!-- -->

# Genotype conditional probabilities

In order to use the genetic map in QTL, we need to obtain the conditional probability of all possible 400 genotypes along the three linkage groups for all individuals in the full-sib population. This can be computed using the function `calc_genoprob`. 


```r
genoprob <- lapply(MAPs, calc_genoprob)
```

Each position of the object `genoprob` contains two elements: an array of dimensions $400 \times number \; of \; markers \times  number \; of \; individuals$ and the position of the markers in the maps in centimorgans. A graphical representation of the genotype probabilities along the three linkage groups in any individual (in this case individual 1) can be obtained using


```r
ind <- 1
dg <- sapply(genoprob, function (x) max(x$map))
dg <- dg/max(dg)
layout(matrix(1:12, ncol = 4), widths = dg)
for(i in 1:12)
{
  d <- genoprob[[i]]$map
  image(t(genoprob[[i]]$probs[,,ind]),
        col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
        axes=FALSE,
        xlab = "Markers",
        ylab = "",
        main = paste("LG", i))
  axis(side = 1, at = d/max(d),
       labels =rep("", length(d)), las=2)
}
```

![](README_files/figure-html/plot_genoprob_bla-1.png)<!-- -->

In this figure, the x-axis represents the genetic map and the y-axis represents the 400 possible genotypes in the full-sib population. The color scale varies from dark purple (high probabilityes) to light yellow (low probabilities). The `genoprob` object obtained here can be used to perform QTL analysis using the R package `QTLpoly` [@Pereira2019], which is an under development software to map multiple QTLs in full-sib families of outcrossing autopolyploid species. 


```r
genoprob.err <- vector("list", 12)
for(i in 1:12)
genoprob.err[[i]] <- calc_genoprob_dist(input.map = MAPs.err[[i]], dat.dist =dat.pot.filt)
```


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

![](README_files/figure-html/plot_genoprob_err-1.png)<!-- -->



# References
