---
title: Building a genetic map using potato genotype data from SolCAP (with minimal comments)
author: "Marcelo Mollinari"
date: "2018-10-23"
output:
 html_document:
   highlight: tango
   keep_md: yes
   theme: united
   toc: yes
   toc_float:
     collapsed: no
 md_document:
   variant: markdown_github
 pdf_document:
   toc: yes
   toc_depth: '3'
   highlight: tango
linestretch: 1.2
#bibliography: biblio.bib
#output: rmarkdown::html_vignette
#vignette: >
#  %\VignetteEngine{knitr::rmarkdown}
---



# Context

The Solanaceae Coordinated Agricultural Project (SolCAP) consitium provides sequences and SNPs used in the SolCAP 8300 Infinium Chip in their [website](http://solcap.msu.edu/potato_infinium.shtml). Here, we will use the first 4x mapping population in Coded Potato Infinium GenomeStudio Project data set. 

# Genotype calling using ClusterCall 

NOTE: Here, we use the marker dosage with no associated probability distribution. 


```r
require(ClusterCall)
require(mappoly)
setwd("~/repos/tutorials/solcap/")
datxy<-read.table(file = "~/repos/tutorials/solcap/potato_x_y_solcap_public.txt", header=TRUE)
datxy[1:10,1:10]
colnames(datxy)
X<-grep(pattern = ".X", colnames(datxy))
Y<-grep(pattern = ".Y", colnames(datxy))
datx<-datxy[,X]
daty<-datxy[,Y]

Px<-apply(datx[,1:2], 1, sum)
Py<-apply(daty[,1:2], 1, sum)
Qx<-apply(datx[,3:4], 1, sum)
Qy<-apply(daty[,3:4], 1, sum)

raw2theta<-function(x,y)
  2*atan(y/x)/pi

raw2r<-function(x,y)
  sqrt(x^2 + y^2)

## Theta
dat.theta<-matrix(NA, nrow(datx), ncol(datx)-2)
dat.theta[,1]<-raw2theta(Px, Py)
dat.theta[,2]<-raw2theta(Qx, Qy)
for(i in 5:ncol(datx))
{
  dat.theta[,i-2]<-raw2theta(datx[,i], daty[,i])
}
dimnames(dat.theta)<-list(as.character(datxy[,1]), c("P1", "P2", unlist(strsplit(colnames(datx), ".X")[-c(1:4)])))
dim(dat.theta)
dat.theta[1:10, 1:10]
write.csv(x = dat.theta, file = "~/repos/tutorials/solcap/solcap_public_theta.csv")

## r
dat.r<-matrix(NA, nrow(datx), ncol(datx)-2)
dat.r[,1]<-raw2r(Px, Py)
dat.r[,2]<-raw2r(Qx, Qy)
for(i in 5:ncol(datx))
{
  dat.r[,i-2]<-raw2r(datx[,i], daty[,i])
}
dim(dat.r)
dimnames(dat.r)<-list(as.character(datxy[,1]), c("P1", "P2", unlist(strsplit(colnames(datx), ".X")[-c(1:4)])))
write.csv(x = dat.r, file = "~/repos/tutorials/solcap/solcap_public_r.csv")

###Running ClusterCall
ab <- read.pop(theta.file = "solcap_public_theta.csv", r.file = "solcap_public_r.csv")
AB <- CC.bipop(ab, parent1 = "P1", parent2 = "P2", n.core = 16)
mrks<-rownames(ab@theta)
i<-1
inspect.marker(AB, mrks[i])
inspect.marker(AB, "solcap_snp_c2_36691")
inspect.marker(AB, "solcap_snp_c2_36686")

abline(h = by(data = AB@theta[mrks[i], ], INDICES = AB@geno[mrks[i], ], mean),
       lty = 2, lwd = 0.5, col = 2)
dev.off()

## Genomic infomation
blast <- read.csv(file="~/repos/potato_analysis_to_publish/blast_solcap.out", header=TRUE)
## selecting the best hits
blast <-blast[!duplicated(blast$query.id),]
dim(blast)
head(blast)
## filter the markers based on the identity ( > 95.5)
blast$snp.pos<-blast$s..start+ceiling((blast$s..start - blast$s..end)/2)
plot(sort(blast$"X..identity"), xlab = "SNP", ylab = "Identity")
abline(h=95.5)
dev.off()
## percentage discarded
100*sum(blast$"X..identity" < 95.5)/nrow(blast)
blast<-blast[blast$"X..identity"> 95.5,]
## Creating a list with all SNPs allocated in the 12 potato
## chromossomes and their positions
chrnames<-sort(as.character(unique(blast$subject.id)))
chr<-vector("list", length(chrnames))
names(chr)<-chrnames
## Attributing chromosomes and positions to SNPs
for(i in chrnames){
  #i<-"chr01"
  temp<-subset(blast, subject.id == i)
  chr[[i]]<-temp[order(temp$snp.pos),c("query.id", "snp.pos")]
}
sum(unlist(sapply(chr, nrow))[1:12])
nrow(blast)
## allowed percentage of double reduction
dr.limit <- 5
## allowed percentage of missing data
miss.limit <- 10
pos<-ch<-nm<-geno<-NULL
for(k in 1:length(chr))
{
  cat("\n Chromosome ", k, " :", sep="")
  for(j in 1:nrow(chr[[k]]))
  {
    i<-as.character(chr[[k]]$query.id[j])
    if(sum(AB@info$marker%in%i) != 1)
      next()
    new.name<-paste(i, names(chr[k]), round(chr[[k]]$snp.pos[j]/1000000,2), sep="_")
    parent.dose<-c(AB@geno[i,AB@parent1], AB@geno[i,AB@parent2])
     if(any(is.na(parent.dose)))
    {
      #print(i)
      next()
    }
    if(all(parent.dose==c(4,0)) | all(parent.dose==c(0,4)) | all(parent.dose==c(0,0)) | all(parent.dose==c(4,4)))
      next()
    #print(parent.dose)
    if(length(table(AB@geno[i,])) ==1)
      next()
    g<-which(segreg_poly(m = 4, AB@geno[i,AB@parent1], AB@geno[i,AB@parent2]) > 0) - 1
    dr<-!AB@geno[i, -c(1:2)]%in%g
    ## if the number os individuals with double reduction
    ## is bigger than the threshold, discard the marker
    if(sum(dr) > ceiling(dr.limit * (ncol(ab@theta) - 2)/100))
      next()
    ## eliminating markers with more missing data than the allowed threshold
    if(sum(is.na(AB@geno[i, -c(1:2)])) >  ceiling(miss.limit * (ncol(ab@theta) - 2)/100))
      next()
    if(sum(dr) != 0)
      AB@geno[i, c(FALSE, FALSE, dr)]<-NA
    #inspect.marker(AB, i)
    #abline(h = by(data = AB@theta[i, ], INDICES = AB@geno[i, ], mean),
    #       lty = 2, lwd = 0.5, col = 2)
    p<-segreg_poly(m = 4, dP = AB@geno[i,1], dQ = AB@geno[i,2])
    v2<-rep(0,5)
    names(v2)<-names(p)<-0:4
    v<-table(AB@geno[i,])
    v2[names(v)]<-v
    test<-chisq.test(v2[p!=0],  p = p[p!=0])
    if(test$p.value < 0.05)
      next()
    nm<-c(nm, new.name)
    #geno.temp<-AB@geno[i,]#-c(1:2)]
    #geno.temp[!geno.temp%in%names(p)[p!=0]]<-NA
    geno<-rbind(geno, AB@geno[i,])
    ch<-c(ch, k)
    pos<-c(pos, chr[[k]]$snp.pos[j])
    #cat(".")
  }
}
rownames(geno)<-nm
dim(geno)
colnames(geno)[-c(1:2)]<-sapply(strsplit(x = colnames(geno), split = "X4x_population1"), function(x) paste0("ind", x[2]))[-c(1:2)]
geno[1:10, 1:10]
## ClusterCall to polymap2
indnames<-colnames(geno)[-c(1:2)]
mrknames<-rownames(geno)
write(paste("ploidy", 4), file="potato_solcap")
write(paste("nind", length(indnames)), file="potato_solcap", append=TRUE)
write(paste("nmrk", nrow(geno)), file="potato_solcap", append=TRUE)
cat("mrknames", mrknames, file="potato_solcap", append=TRUE)
cat("\nindnames", indnames, file="potato_solcap", append=TRUE)
cat("\ndosageP", geno[,AB@parent1], file="potato_solcap", append=TRUE)
cat("\ndosageQ", geno[,AB@parent2], file="potato_solcap", append=TRUE)
cat("\nseq", ch, file="potato_solcap", append=TRUE)
cat("\nseqpos", pos, file="potato_solcap", append=TRUE)
write("\nnphen 0", file="potato_solcap", append=TRUE)
write("pheno---------------------------------------", file="potato_solcap", append=TRUE)
write("geno---------------------------------------", file="potato_solcap", append=TRUE)
dose_to_dist<-function(x, m=4, p, q)
{
  if(is.na(x))
    return(segreg_poly(m, dP = p, dQ = q))
  y<-numeric(m+1)
  y[x+1]<-1
  return(y)
}
geno.out<-NULL
for(i in 1:length(indnames))
{
  a<-t(mapply(dose_to_dist, geno[,i+2], m=4, p=geno[,AB@parent1], q=geno[,AB@parent2]))
  colnames(a)<-NULL
  geno.out<-rbind(geno.out, data.frame(rownames(a), indnames[i], a))
}
write.table(geno.out, file="potato_solcap", append=TRUE, quote=FALSE,
            row.names=FALSE, col.names=FALSE)
```

# Map construction 

## Reading data set


```r
require(mappoly)
```

```
## Loading required package: mappoly
```

```r
solcap.dat<-read_geno_dist(file.in = "~/repos/tutorials/solcap/potato_solcap")
```

```
## Reading the following data:
##     Ploidy level: 4
##     No. markers:  3712
##     No. individuals:  160
##     This dataset contains sequence information.
##    ...
##     Done with reading.
```

```r
solcap.dat
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:    4 
##     No. individuals:    160 
##     No. markers:        3712 
## 
##     This dataset contains sequence information.
##     ----------
##     No. of markers per dosage in both parents:
##      P Q freq
##      0 1  280
##      0 2  124
##      0 3   37
##      1 0  284
##      1 1  229
##      1 2  184
##      1 3  101
##      1 4   23
##      2 0   93
##      2 1  167
##      2 2  223
##      2 3  180
##      2 4  144
##      3 0   39
##      3 1  119
##      3 2  241
##      3 3  320
##      3 4  361
##      4 1   34
##      4 2  152
##      4 3  377
```

Summary 


```r
print(solcap.dat, detailed = TRUE)
```

```
## This is an object of class 'mappoly.data'
##     Ploidy level:    4 
##     No. individuals:    160 
##     No. markers:        3712 
## 
##     ----------
##     No. markers per sequence:
##        seq No.mrk
##          1    400
##          2    284
##          3    238
##          4    426
##          5    284
##          6    425
##          7    410
##          8    297
##          9    265
##         10    173
##         11    289
##         12    221
##     ----------
##     Markers not allocated: 0
##     ----------
##     No. of markers per dosage in both parents:
##      P Q freq
##      0 1  280
##      0 2  124
##      0 3   37
##      1 0  284
##      1 1  229
##      1 2  184
##      1 3  101
##      1 4   23
##      2 0   93
##      2 1  167
##      2 2  223
##      2 3  180
##      2 4  144
##      3 0   39
##      3 1  119
##      3 2  241
##      3 3  320
##      3 4  361
##      4 1   34
##      4 2  152
##      4 3  377
```

```r
plot(solcap.dat)
```

![](solcap_example_files/figure-html/dataset-1.png)<!-- -->

### Recombination fraction estimation

Now, let us compute the recombination fraction between all markers in the data set. 


```r
s<-make_seq_mappoly(solcap.dat, arg = "all")
counts.web<-cache_counts_twopt(s, get.from.web = TRUE)
```

```
## Internet conectivety ok.
## Loading genotype counts from web
```

```r
all.pairs<-est_pairwise_rf(input.seq = s,
                           count.cache = counts.web,
                           n.clusters = 16, #toke approx. 9 minutes using 24 CPUs in BRC cluster
                           verbose=TRUE) 
```

```
## INFO: Using  16  CPUs for calculation.
## INFO: Done with 6887616  pairs of markers 
## INFO: Calculation took: 597.696 seconds
```

```r
mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
```

```
## INFO: Going singlemode. Using one CPU.
```

```r
mat.full
```

```
##   This is an object of class 'mappoly.rf.matrix'
## 
##   Criteria used to filter markers:
## 
##       Configuration phase LOD:            0 
##       Recombination fraction LOD:         0 
##       Maximum recombination fraction:     0.5 
## 
##   No. markers:              3712 
##   Percentage filled:        86.2 %
```

```r
plot(mat.full)
```

![](solcap_example_files/figure-html/rf_estimation-1.png)<!-- -->

### Grouping


```r
lgs<-group_mappoly(input.mat = mat.full,
                   input.seq = s,
                   expected.groups = 13,
                   comp.mat = TRUE, 
                   inter = FALSE)
print(lgs, detailed = TRUE)
```

```
##   This is an object of class 'mappoly.group'
##   ------------------------------------------
##   Criteria used to assign markers to groups:
## 
##     - Number of markers =         3712 
##     - Number of linkage groups = 12 
##     - Number of markers per linkage groups: 
##     group n.mrk
##         1   418
##         2   286
##         3   252
##         4   432
##         5   285
##         6   415
##         7   410
##         8   291
##         9   270
##        10   162
##        11   284
##        12   207
##   ------------------------------------------
##      1   2   3   4   5   6   7   8   9  10  11  12 NH
## 1  390   0   9   5   2   1   0   2   1   0   3   5  0
## 2    2 202   0   0   0   0   0   0   0   3   0   0  0
## 3    1   3 268   0   3   1   0   0   0   5   3   1  0
## 4    1   1   0 267   0   4   7   2   1   0   2   0  0
## 5    1   7   3   0 229   1   1   5   0   0   1   4  0
## 6    2   0   0   3   2 395   0   2   0   6   0   0  0
## 7    1   4   0   1   0   1 254   1   2   0   5   1  0
## 8    1   0   0   1   0   5   1 154   0   0   0   0  0
## 9    1   1   3   1   1   1   1   0 420   1   1   1  0
## 10   0   0   1   0   0   0   0   0   0   0   0   0  0
## 11   0   0   0   1   1   0   1   2   2 281   2   1  0
## 12   0   1   0   3   0   1   0   5   0   0 404   1  0
## 13   0   2   0   2   0   0   0   0   0   1   4 275  0
## 
##   ------------------------------------------
```

```r
plot(lgs)
```

![](solcap_example_files/figure-html/grouping-1.png)<!-- -->

```r
## Selecting groups
LGS<-lapply(c(1:12), function(x, lgs) make_seq_mappoly(lgs, x), lgs)
```

### Filtering and ordering markers within groups


```
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
## INFO: Going singlemode. Using one CPU.
```

```
## [[1]]
## [1] TRUE
## 
## [[2]]
## [1] TRUE
## 
## [[3]]
## [1] TRUE
## 
## [[4]]
## [1] TRUE
## 
## [[5]]
## [1] TRUE
## 
## [[6]]
## [1] TRUE
## 
## [[7]]
## [1] TRUE
## 
## [[8]]
## [1] TRUE
## 
## [[9]]
## [1] TRUE
## 
## [[10]]
## [1] TRUE
## 
## [[11]]
## [1] TRUE
## 
## [[12]]
## [1] TRUE
```

Now, let us plot the reordered recombination fraction 

```r
op <- par(mfrow = c(3, 4), pty = "s")
for(i in 1:12)
  plot(M.filt[[i]], ord = MDS.seq[[i]]$seq.mrk.names, 
       main.text = paste0("LG", i), index = FALSE)
```

![](solcap_example_files/figure-html/mds_plot-1.png)<!-- -->

```r
par(op)
```

Now, the full recombination fraction matrix

```r
o<-unlist(sapply(MDS.seq, function(x) x$seq.mrk.names))
plot(mat.full, ord = o)
```

![](solcap_example_files/figure-html/full_rf_mat-1.png)<!-- -->

### Phasing and multilocus reconstruction

Given the MDS order, we estimate the phase and the recombination fraction between all adjacent markers. After that, we reestimate the map using the function `est_full_hmm_with_global_error`, which includes a global genotype error in the HMM model.


```r
## function for parallel computation
my.phase.func<-function(X)
{
  return(est_rf_hmm_sequential(input.seq = X[[1]],
                               thres.twopt = 10,
                               thres.hmm = 10,
                               extend.tail = 50,
                               twopt = X[[2]],
                               verbose = TRUE,
                               tol = 10e-3,
                               tol.final = 10e-4,
                               phase.number.limit = 40,
                               sub.map.size.diff.limit = 6,
                               info.tail = TRUE,
                               reestimate.single.ph.configuration = TRUE, 
                               high.prec = FALSE))
  
}
## list with ordered sequences and pairwise recombination fractions
X<-vector("list", 12)
for(i in 1:12)
  X[[i]]<-list(MDS.seq[[i]], P.filt[[i]])

## Running phasing algorithm in 12 CPUs
system.time({
  cl <- parallel::makeCluster(12)
  parallel::clusterEvalQ(cl, require(mappoly))
  parallel::clusterExport(cl, "solcap.dat")
  MAPS <- parallel::parLapply(cl,X,my.phase.func)
  MAPs.res <- parallel::parLapply(cl, MAPS, est_full_hmm_with_global_error, 
                                  error = 0.05, tol = 10e-4)
  parallel::stopCluster(cl)
})
```

```
##     user   system  elapsed 
##    2.904    0.637 1579.389
```

```r
MAPs.res
```

```
## [[1]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 362 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 155.09
##        log-likelihood:	 -9682.29
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[2]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 249 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 119.95
##        log-likelihood:	 -7641.79
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[3]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 223 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 100.46
##        log-likelihood:	 -6728.69
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[4]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 377 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 117.98
##        log-likelihood:	 -9944.5
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[5]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 250 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 93.01
##        log-likelihood:	 -7182.98
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[6]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 357 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 95.13
##        log-likelihood:	 -7983.67
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[7]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 365 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 102.81
##        log-likelihood:	 -8321.53
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[8]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 254 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 90.54
##        log-likelihood:	 -6500.04
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[9]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 238 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 97.18
##        log-likelihood:	 -6763.53
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[10]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 135 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 107.39
##        log-likelihood:	 -5152.62
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[11]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 253 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 98.39
##        log-likelihood:	 -6914.15
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
## 
## [[12]]
## This is an object of class 'mappoly.map'
##     Ploidy level:	 4 
##     No. individuals:	 160 
##     No. markers:	 177 
##     No. linkage phases:	 1 
## 
##     ---------------------------------------------
##     Number of linkage phase configurations:  1
##     ---------------------------------------------
##     Linkage phase configuration:  1
##        map length:	 89.55
##        log-likelihood:	 -5590.01
##        LOD:		 0
##     ~~~~~~~~~~~~~~~~~~
```

```r
sapply(MAPs.res, plot)
```

![](solcap_example_files/figure-html/phasing_and_reestimate-1.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-2.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-3.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-4.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-5.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-6.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-7.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-8.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-9.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-10.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-11.png)<!-- -->![](solcap_example_files/figure-html/phasing_and_reestimate-12.png)<!-- -->

```
##      [,1]        [,2]        [,3]        [,4]         [,5]        
## path NULL        NULL        NULL        NULL         NULL        
## name "GRID.VP.3" "GRID.VP.6" "GRID.VP.9" "GRID.VP.12" "GRID.VP.15"
## n    1           1           1           1            1           
##      [,6]         [,7]         [,8]         [,9]         [,10]       
## path NULL         NULL         NULL         NULL         NULL        
## name "GRID.VP.18" "GRID.VP.21" "GRID.VP.24" "GRID.VP.27" "GRID.VP.30"
## n    1            1            1            1            1           
##      [,11]        [,12]       
## path NULL         NULL        
## name "GRID.VP.33" "GRID.VP.36"
## n    1            1
```




