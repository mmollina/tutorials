require(mappoly)
require(polymapR)
source("~/repos/MAPPoly/docs/polymapR_example/polymapR_to_mappoly.R")
data("segregating_data")
tetra.data<-polymapR_to_mappoly(segregating_data)
print(tetra.data, detailed = TRUE)

s<-make_seq_mappoly(tetra.data, arg = "all")
counts.web<-cache_counts_twopt(s, get.from.web = TRUE)
all.pairs<-est_pairwise_rf(input.seq = s,
                           count.cache = counts.web,
                           n.clusters = 8,
                           verbose=TRUE)
mat.full<-rf_list_to_matrix(input.twopt=all.pairs, n.clusters = 5)
lgs<-group_mappoly(input.mat = mat.full,
                   input.seq = s,
                   expected.groups = 12,
                   comp.mat = FALSE, 
                   inter = TRUE)
lgs
LGS<-lapply(c(1,2,4,5,6), function(x, lgs) make_seq_mappoly(lgs, x), lgs)
P<-lapply(LGS, make_pairs_mappoly, input.twopt = all.pairs)
LGS.filt<-lapply(P, rf_snp_filter, thresh.LOD.ph = 5, 
                 thresh.LOD.rf = 5, thresh.rf = 0.15, 
                 thresh.perc = 0.10)
P.filt<-lapply(LGS.filt, make_pairs_mappoly, input.twopt = all.pairs)
M.filt<-lapply(P.filt, rf_list_to_matrix)
MDS.ord<-lapply(M.filt, mds_mappoly)
MDS.seq<-lapply(MDS.ord, make_seq_mappoly)

op <- par(mfrow = c(2, 3), pty = "s")
for(i in 1:5)
  plot(M.filt[[i]], ord = MDS.seq[[i]]$seq.mrk.names, 
       main.text = paste0("LG", i), index = FALSE)
par(op)

o<-unlist(sapply(MDS.seq, function(x) x$seq.mrk.names))
plot(mat.full, ord = o)

load("~/repos/MAPPoly/docs/polymapR_exemple.RData")

system.time(
  {
    id<-c(4,3,1,2,5)
    pal<-comp.PHQ<-comp.PHP<-MAPS<-vector("list", 5)
    for(i in 2:5)
    {
      MAPS[[i]]<-est_rf_hmm_sequential(MDS.seq[[i]],
                                       thres.twopt = 10,
                                       thres.hmm = 10,
                                       extend.tail = 20,
                                       twopt = P[[i]],
                                       verbose = TRUE,
                                       tol = 10e-3,
                                       tol.final = 10e-5,
                                       phase.number.limit = 20,
                                       sub.map.size.diff.limit = 4,
                                       info.tail = TRUE,
                                       reestimate.single.ph.configuration = TRUE, 
                                       high.prec = FALSE)
      
      YP<-phased.maplist[[id[i]]][,3:6]
      rownames(YP)<-phased.maplist[[id[i]]][,1]
      YQ<-phased.maplist[[id[i]]][,7:10]
      rownames(YQ)<-phased.maplist[[id[i]]][,1]
      mrk.names<-tetra.data$mrk.names[MAPS[[i]]$maps[[1]]$seq.num]
      nm<-intersect(rownames(YQ), mrk.names)
      AP<-MAPS[[i]]$maps[[1]]$seq.ph$P
      AQ<-MAPS[[i]]$maps[[1]]$seq.ph$Q
      names(AQ)<-names(AP)<-mrk.names
      BP<-ph_matrix_to_list(YP)
      names(BP)<-rownames(YP)
      BQ<-ph_matrix_to_list(YQ)
      names(BQ)<-rownames(YQ)
      comp.PHP[[i]]<-compare_haplotypes(m = 4, h1 = AP[nm], h2 = BP[nm]) 
      comp.PHQ[[i]]<-compare_haplotypes(m = 4, h1 = AQ[nm], h2 = BQ[nm]) 
      pal[[i]]<-mrk.names%in%rownames(YQ)*2+2
      MAPs.res[[i]]<-est_full_hmm_with_global_error(input.map = MAPS[[i]], error = 0.05, tol = 10e-4)
    }
  }
)

a<-phased.maplist$LG4[,1:2]
b<-data.frame(marker = tetra.data$mrk.names[MAPS[[i]]$maps[[1]]$seq.num],
              position = cumsum(imf_h(c(0, MAPS[[i]]$maps[[1]]$seq.rf))))
head(a);head(b)
plot(match(a$marker, b$marker))


s.comp1<-make_seq_mappoly(tetra.data, arg = intersect(a$marker, b$marker))
s.comp2<-make_seq_mappoly(tetra.data, arg = intersect(b$marker, a$marker))
match(s.comp1$seq.num, s.comp2$seq.num)


oa<-match(s.comp1$seq.num, MAPs.res$maps[[1]]$seq.num)
ob<-match(s.comp2$seq.num, MAPs.res$maps[[1]]$seq.num)
ares<-get_submap(MAPs.res[[1]], mrk.pos = oa, reestimate.rf = TRUE)
bres<-get_submap(MAPs.res[[1]], mrk.pos = ob, reestimate.rf = TRUE)
ares; bres

areserr<-est_full_hmm_with_global_error(ares, error = 0.05, tol = 10e-4)
breserr<-est_full_hmm_with_global_error(bres, error = 0.05, tol = 10e-4)
areserr; breserr

plot(MAPs.res[[1]], col.cte = pal[[1]])





 