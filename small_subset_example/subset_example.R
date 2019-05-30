require(mappoly)
data("tetra.solcap")

## selecting SNPs in the first chromosome
seq1<-make_seq_mappoly(tetra.solcap, arg = 'seq1')


## subset based on marker positions
#####
set.seed(1234)
(subset.1<-sample(seq1$seq.num, 100)) ## sampling 100 markers
(subset.1<-subset.1[order(seq1$sequence.pos[subset.1])]) ## ordering according the genome info
seq.subset.1<-make_seq_mappoly(input.obj = tetra.solcap, arg = subset.1) ## making a sequence
plot(seq.subset.1)

## subset based on marker names
#####
(subset.2<-tetra.solcap$mrk.names[subset.1]) ## sampling 100 markers
seq.subset.2<-make_seq_mappoly(input.obj = tetra.solcap, arg = subset.1) ## making a sequence
plot(seq.subset.2)

identical(seq.subset.1, seq.subset.2)

## physical vs. mrk order
plot(seq.subset.1$seq.num~seq.subset.1$sequence.pos) 
## Two-point recombination fraction estimation
counts<-cache_counts_twopt(seq.subset.1, 
                           get.from.web=TRUE)
subset.pairs<-est_pairwise_rf(seq.subset.1,
                              count.cache = counts,
                              n.clusters=1)
plot(rf_list_to_matrix(subset.pairs, thresh.LOD.ph = 3, thresh.LOD.rf = 3), index = FALSE)
## Map construction
subset.map <- est_rf_hmm_sequential(input.seq  = seq.subset.1,
                                    start.set = 5,
                                    thres.twopt = 10, 
                                    thres.hmm = 10,
                                    extend.tail = 20,
                                    info.tail = TRUE, 
                                    twopt = subset.pairs,
                                    phase.number.limit = 20,
                                    sub.map.size.diff.limit = 10,
                                    reestimate.single.ph.configuration = TRUE,
                                    tol = 10e-3,
                                    tol.final = 10e-4)
subset.map
plot(subset.map)

## Reestimating map with global error = 0.05
subset.map.reest<-est_full_hmm_with_global_error(input.map = subset.map, 
                                                 error=0.05, 
                                                 tol=10e-4, 
                                                 verbose = TRUE)
subset.map.reest
plot(subset.map.reest)

## Compiting posterior state probabilities
op<-par(mfrow = c(1,2))
genoprob.lg1<-calc_genoprob(input.map = subset.map) ## no global error
genoprob.lg1.prior<-calc_genoprob(input.map = subset.map.reest) ## considering global error
ind <- 1
d <- genoprob.lg1$map
pr<-genoprob.lg1$probs[,,ind]
image(t(pr),
      col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
      axes=FALSE,
      xlab = "Markers",
      ylab = "",
      main = paste("LG_1 - ind ", ind))
axis(side = 1, at = d/max(d),
     labels =rep("", length(d)), las=2)
axis(side = 2, at = seq(0,1,length.out = nrow(pr)),
     labels = rownames(pr), las=2, cex.axis=.5)

d <- genoprob.lg1.prior$map
pr<-genoprob.lg1.prior$probs[,,ind]
image(t(pr),
      col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
      axes=FALSE,
      xlab = "Markers",
      ylab = " ",
      main = paste("LG_1 - ind ", ind, "w/ error"))
axis(side = 1, at = d/max(d),
     labels =rep("", length(d)), las=2)



