## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache = TRUE)
knitr::opts_chunk$set(eval = FALSE)


## ----load_mappoly, eval=TRUE, results='hide'-----------------------------
library(mappoly)


## ---- data set_example, out.width='40%', fig.align='center', fig.cap='Figure 1: Example of data set containing genotype probabilities', eval=TRUE, echo=FALSE----
knitr::include_graphics('image/dataset_tetra.png')


## ----load_data, eval=TRUE------------------------------------------------
#(~46 seconds)
dat.pot<-read_geno_dist(file.in = "SolCAP", prob.thres = 0.9)


## ----print_datset, eval=TRUE---------------------------------------------
print(dat.pot, detailed = TRUE)


## ----plot_datset, eval=TRUE----------------------------------------------
plot(dat.pot)


## ---- filtering_missing, eval=TRUE---------------------------------------
dat.pot.filt<-filter_missing(dat.pot, filter.thres = 0.2, inter = FALSE)
print(dat.pot.filt)
plot(dat.pot.filt) #based in a bonferroni aproximation


## ---- select_mrk_chisq, eval = TRUE--------------------------------------
seq.filt<-filter_segregation(dat.pot.filt, chisq.pval.thres = 0.05/dat.pot.filt$n.mrk, inter = FALSE)
seq.filt<-make_seq_mappoly(seq.filt)
print(seq.filt)


## ---- plot_mrk_info, eval = TRUE, out.width = "1000px"-------------------
plot_mrk_info(input.data = dat.pot.filt, mrk = "solcap_snp_c1_13686")


## ---- redundant_elimi, eval=TRUE-----------------------------------------
## Filtering out redundant markers
seq.uni<-elim_redundant(input.seq = seq.filt)
plot(seq.uni)


## ---- select_all, eval=TRUE----------------------------------------------
## Filtering out redundant markers
seq.all.lgs<-make_seq_mappoly(seq.uni)


## ---- load_counts, eval=TRUE---------------------------------------------
counts<-cache_counts_twopt(input.seq = seq.all.lgs, get.from.web = TRUE)
counts


## ---- cache_twopt, eval=TRUE, echo=FALSE---------------------------------
load("all_rf_pairwise_for_3654_mrks.RData")


## ----two_pt, eval=FALSE--------------------------------------------------
## #(~ 9.5 minutes)
## all.rf.pairwise <- est_pairwise_rf(input.seq = seq.all.lgs,
##                                    count.cache = counts,
##                                    n.clusters = 16)
## all.rf.pairwise


## ---- plot_twopt_example, eval=TRUE, out.width = "500px", fig.align="center"----
all.rf.pairwise$pairwise$`93-98`
plot(all.rf.pairwise, first.mrk = 93, second.mrk = 98)


## ----rf_mat, echo=TRUE, eval = TRUE--------------------------------------
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)


## ----plot_full_mat, dpi=200, eval=TRUE-----------------------------------
id<-get_genomic_order(seq.all.lgs)
plot(mat, ord = rownames(id), index = FALSE)


## ----group, eval=TRUE----------------------------------------------------
grs <- group_mappoly(input.mat = mat,
                     input.seq = seq.all.lgs,
                     expected.groups = 12,
                     comp.mat = TRUE, 
                     inter = TRUE)
grs


## ----plot_group, eval = TRUE---------------------------------------------
plot(grs)


## ----make_lgs, eval=TRUE-------------------------------------------------
LGS<-vector("list", 12)
for(j in 1:12){
  temp1<-make_seq_mappoly(grs, j)
  temp2<-get_genomic_order(temp1) # assembling sequence considering the genomic order
  nm<-names(which(temp2[,1]==names(which.max(table(temp2[,1])))))
  lgtemp<-make_seq_mappoly(dat.pot.filt, nm)
  LGS[[j]]<-list(lg = lgtemp, 
                 tpt = make_pairs_mappoly(all.rf.pairwise, input.seq = lgtemp))
}


## ---- all_mat_rf, eval=TRUE, echo = FALSE, results='hide'----------------
op <- par(mfrow = c(3, 4), pty = "s", mar=c(1,1,1,1))
for(i in 1:12)
  plot(rf_list_to_matrix(LGS[[i]]$tpt), ord = LGS[[i]]$lg$seq.mrk.names, 
       main.text = paste0("LG", i), index = FALSE)
par(op)


## ---- map_lg12, eval=TRUE, results='hide'--------------------------------
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


## ---- map_lg12_plot, echo=TRUE, eval=TRUE--------------------------------
print(lg12.map)
plot(lg12.map)


## ---- map_err_lg12, eval=TRUE, results='hide'----------------------------
lg12.map.prior<-reest_map_geno_dist(lg12.map, dat.dist = dat.pot.filt, verbose = FALSE, tol = 10e-4)
plot(lg12.map.prior$reestimated.map)


## ---- map_prior_lg12, eval=TRUE------------------------------------------
lg12.map.error<-est_full_hmm_with_global_error(input.map = lg12.map, error = 0.01)
plot(lg12.map.error)


## ---- genoprob_12_hard, eval = TRUE--------------------------------------
genoprob.lg12<-calc_genoprob(input.map = lg12.map)


## ----plot_genoprob, eval=TRUE,  out.width='200%',,fig.align="center"-----
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


## ---- genoprob_12_soft, eval = TRUE--------------------------------------
genoprob.lg12.prior<-calc_genoprob_dist(input.map = lg12.map.error, dat.dist = dat.pot.filt)


## ----plot_genoprob_prior, eval=TRUE--------------------------------------
ind <- 1
d <- genoprob.lg12.prior$map
pr<-genoprob.lg12.prior$probs[,,ind]
image(t(pr),
      col=RColorBrewer::brewer.pal(n=9 , name = "YlOrRd"),
      axes=FALSE,
      xlab = "Markers",
      ylab = " ",
      main = paste("LG", i))
axis(side = 1, at = d/max(d),
     labels =rep("", length(d)), las=2)
axis(side = 2, at = seq(0,1,length.out = nrow(pr)),
     labels = rownames(pr), las=2, cex.axis=.5)


## ----mat_for_ordering, results='hide', eval=TRUE-------------------------
mt <- lapply(LGS, function(x) rf_list_to_matrix(x$tpt))


## ----mds, results='hide', eval=TRUE--------------------------------------
mds.ord <- lapply(mt, mds_mappoly)


## ----compare_order, eval=TRUE, results='hide'----------------------------
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
par(op)


## ---- map_lg12_mds, eval=TRUE, results='hide'----------------------------
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


## ----lg12_map_mds_plot, results='hide', eval=TRUE------------------------
plot(lg12.map.mds)
print(lg12.map)
print(lg12.map.mds)


## ----map_comp, eval=TRUE, results='hide'---------------------------------
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


## ----hmm_map, results='hide', eval = TRUE--------------------------------
## Performing parallel computation
#(~20 minutes)
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


## ----all_maps_error, results='hide', eval = TRUE-------------------------
MAPs.err <- MAPs.prior <- vector("list", 12)
for(i in 1:12){
  cat("LG:", i ,"\n")
   MAPs.err[[i]] <- est_full_hmm_with_global_error(input.map = MAPs[[i]], error = 0.05, tol = 10e-4, verbose = FALSE)
   MAPs.prior[[i]] <- reest_map_geno_dist(input.map = MAPs[[i]], dat.dist = dat.pot.filt, verbose = FALSE, tol = 10e-4)
}


## ----print_maps, results='hide'------------------------------------------
sapply(MAPs, function(x) plot(x))
sapply(MAPs.err, function(x) plot(x))
sapply(MAPs.prior, function(x) plot(x$reestimated.map))


## ----genoprob_2, eval=TRUE, results='hide'-------------------------------
genoprob.err <- vector("list", 12)
for(i in 1:12)
genoprob.err[[i]] <- calc_genoprob_dist(input.map = MAPs.err[[i]], dat.dist = dat.pot.filt)


## ----plot_genoprob_all, eval = TRUE--------------------------------------
ind <- 2
op <- par(mfrow = c(3, 4), mar=c(1,1,1,1))
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
par(op)


## ------------------------------------------------------------------------
require(mappoly)
genoprob.ls <- genoprob.err
saveRDS(genoprob.ls, file = "genoprob_list.rds")
final.map<-MAPs.err
saveRDS(final.map, file = "final.map.rds")
saveRDS(dat.pot.filt, file = "input.data.rds")


## ----save, eval=FALSE, echo=FALSE, include=FALSE, eval=TRUE--------------
save.image(file = "all_data_after_analysis.rda", compress = TRUE)

