## Declaring functions
homo_prob_tetra<-function(ind, genoprob, parent){
  ## 400 genotype names
  geno.names<-dimnames(genoprob$probs)[[1]]
  ## conditional probability matrix
  pr.mat<-genoprob$probs[,,ind]
  ## choosing parent
  if(parent=="P1"){
    lt<-letters[1:4]
  } else {
    lt<-letters[5:8]
  }
  #####
  ## probability matrix for the homologous chromosomes
  pr.hom<-NULL
  for(i in 1:4)
  {
    y<-apply(pr.mat[grep(pattern = lt[i], geno.names),], 2, sum)
    pr.hom<-rbind(pr.hom, y)
  }
  rownames(pr.hom)<-lt
  return(pr.hom)  
}
plot_homo_prob_tetra<-function(pr.hom, map, ind, parent, title.plot = NULL){
  if(is.null(title.plot))
    title.plot <- parent
  ## Transforming probability matrix for the homologous chromosomes into a dataframe
  lt<-rownames(pr.hom)
  if(grepl(pattern = "a", lt[1]))
    pal<-c("#800000", "#9A6324", "#808000", "#e6194B")
  else if(grepl(pattern = "e", lt[1]))
    pal<-c("#3cb44b", "#469990", "#42d4f4", "#4363d8")
  df.pr<-NULL
  for(i in 1:4)
    df.pr<-rbind(df.pr, data.frame(pr = pr.hom[i,], 
                                   hom = lt[i], 
                                   position = map[colnames(pr.hom)]))
  ## ggplot graphic
  p <- ggplot(df.pr, aes(x = position, y = pr, fill = hom, color  = hom)) +
    geom_density(stat = "identity", alpha = 0.7) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
    ggtitle(title.plot) + facet_grid(rows = vars(hom)) + theme_minimal() + ylab(label = "Homologous probabilty")
  return(p)
}
elim_low_prob_homo_segments_tetra<-function(pr.hom, hom.prob.thresh){
  out<-apply(pr.hom, 2, function(x) {
    y<-rep(0,4)
    y[x > hom.prob.thresh]<-1
    y})
  rownames(out)<-rownames(pr.hom)
  # At least one homologous has pr > hom.prob.thresh
  out <- out[,(apply(out, 2, sum) >= 1)]
  return(out)}
get_segment_length_tetra<-function(pr.hom, map, is.gap = FALSE){
  lt<-rownames(pr.hom)
  A <- pr.hom
  if(!is.gap){
    A[A == 1] <- 2
    A[A == 0] <- 1
    A[A == 2] <- 0
  }
  # if(is.gap){
  #   A<-apply(pr.hom, 2, function(x) {
  #     y<-rep(1,6)
  #     y[x > 0]<-0
  #     y})
  # } else {
  #   A<-apply(pr.hom, 2, function(x) {
  #     y<-rep(0,6)
  #     y[x > 0]<-1
  #     y})
  # }
  dimnames(A)[[1]]<-lt
  B0<-apply(A, 1, function(x) abs(diff(x)))
  D0<-apply(B0, 1, function(x) which(x==1))
  D0 <- B0[which(sapply(D0, length) > 0), ,drop=FALSE]
  B1<-vector("list", 4)
  names(B1)<-lt
  for(j in 1:4){
    x <- c(1, which(diff(A[j,])!=0), ncol(A))
    id<-A[j,x[-length(x)]+1]==0
    x <- cbind(x[-length(x)], x[-1])
    xlen <- NULL
    for(i in 1:nrow(x))
      xlen <- c(xlen, map[x[i,2]] - map[x[i,1]])
    B1[[j]] <- cbind(x, xlen)[id,]
    B1[[j]]<-matrix(B1[[j]], ncol = 3, byrow = FALSE)
    colnames(B1[[j]])<-c("start", "end", "length")
  }
  n.breaks<-matrix(apply(B0, 2, sum), nrow = 1, dimnames = list("ind", lt))
  list(segment.lengths = B1, n.breaks)
}
elim_small_and_incomplete_segments_tetra<-function(pr.hom, segment.lengths, seg.length.thresh){
  lt<-rownames(pr.hom)
  for(j in 1:4){
    for(i in 1:nrow(segment.lengths$segment.lengths[[j]]))
      if(nrow(segment.lengths$segment.lengths[[j]]) > 0)
        if(segment.lengths$segment.lengths[[j]][i,3] < seg.length.thresh)
          pr.hom[,segment.lengths$segment.lengths[[j]][i,1]:segment.lengths$segment.lengths[[j]][i,2]]<-NA
  }
  idrem<-which(apply(pr.hom, 2, function(x) all(is.na(x))))
  if(length(idrem) != 0){
    pr.hom.NA <- pr.hom
    pr.hom.<- pr.hom[ ,-idrem]
  } 
  mrk.names<-pr.hom.res<-NULL
  for(i in 1:ncol(pr.hom)){
    if(!all(is.na(pr.hom[,i])))
      # All three homologous have to have high probabilities, otherwise, we estimate using HMM + high quality markers
      #if(sum(pr.hom[,i], na.rm = TRUE) >= 2 ){
      if(sum(pr.hom[,i], na.rm = TRUE) == 2){
        pr.hom.res <- cbind(pr.hom.res, pr.hom[,i])
        mrk.names<-c(mrk.names, colnames(pr.hom)[i])
      }
  }
  colnames(pr.hom.res)<-mrk.names
  return(pr.hom.res)
}
detect_co_tetra<-function(pr.hom, map, dist.thresh = 5){
  if(is.null(pr.hom))
    return(list(result = "Inconclusive", result.par = "Inconclusive"))
  lt<-rownames(pr.hom)
  ## Detecting crossing overs
  A<-apply(pr.hom, 2, function(x) {
    y<-rep(0,4)
    y[order(x)[1:2]]<-1
    y})
  dimnames(A)[[1]]<-lt
  B<-apply(A, 1, function(x) abs(diff(x)))
  D<-apply(B, 1, function(x) which(x==1))
  D <- B[which(sapply(D, length) > 0), ,drop=FALSE]
  result<-"Zero"
  if(nrow(D)!=0){
    hom.pair <- apply(D, 1, function(x) paste0(lt[which(x==1)], collapse = ""))
    position <- round(map[names(hom.pair)],1) 
    z1 <- apply(D, 1, function(x) {y <- which(x==1); if(length(y)==2) return(y) else return(c(NA, NA))})
    from <- z1[1,]
    to <- z1[2,]
    point.id<-match(names(hom.pair), names(map))
    result<-data.frame(hom.pair, position, from, to, point.id)
  }
  if(as.character(result)[1]=="Zero")
    return(list(result = "Zero", result.par = "Zero"))
  id<-is.na(result$from)
  rem = NULL
  if(any(id)){
    rem<-result[which(id),]
    result<-result[-which(id),]
  }
  if(nrow(result)==0)
    return(list(result = "Inconclusive", result.par = "Inconclusive"))
  rem <- which(diff(result$position) < dist.thresh)
  if(length(rem) > 0){
    for(i in 1:length(rem))
      rem <- c(rem, rem[i] + 1)
    result <- result[-rem,,drop = FALSE]
    if(nrow(result)==0)
      return(list(result = "Inconclusive", result.par = "Inconclusive"))
  }
  return(list(result = result, result.par = result[,1:2]))
}
count_cross_tetra<-function(pr.hom, 
                       map,
                       individual,
                       parent,
                       hom.prob.thresh = 0.80,
                       seg.length.thresh = 10,
                       perc.info = 20,
                       dist.thresh = 2,
                       map.mappoly,
                       dat.mappoly,
                       title.plot = NULL)
{
  plot.1<-plot_homo_prob_tetra(pr.hom = pr.hom, map = map, ind = individual, parent, title.plot = title.plot)
  pr.hom.2<-elim_low_prob_homo_segments_tetra(pr.hom = pr.hom, hom.prob.thresh = hom.prob.thresh)
  plot.2<-plot_homo_prob_tetra(pr.hom = pr.hom.2, map = map, ind = individual, parent, title.plot = title.plot)
  segment.lengths<-get_segment_length_tetra(pr.hom = pr.hom.2, map = map)
  pr.hom.3<-elim_small_and_incomplete_segments_tetra(pr.hom = pr.hom.2, 
                                               segment.lengths = segment.lengths, 
                                               seg.length.thresh = seg.length.thresh)
  if(is.null(pr.hom.3)){
      pr.hom[]<-0
      return(list(p = "Inconclusive", 
                  plot.1 = plot.1, 
                  plot.final = NULL, 
                  n.breaks=segment.lengths[[2]],
                  n.breaks.new=NULL,
                  res.par = "Inconclusive"))
  }
  gap.lengths <- get_segment_length_tetra(pr.hom = pr.hom.3, map = map, is.gap = TRUE)
  pr.hom.4<-elim_small_and_incomplete_segments_tetra(pr.hom = pr.hom.3, 
                                               segment.lengths = gap.lengths, 
                                               seg.length.thresh = seg.length.thresh)
  if(is.null(pr.hom.4)){
    pr.hom[]<-0
    return(list(p = "Inconclusive", 
                plot.1 = plot.1, 
                plot.final = plot_homo_prob_tetra(pr.hom = pr.hom, map = map, ind = individual, parent, title.plot = title.plot), 
                n.breaks=segment.lengths[[2]],
                n.breaks.new=NULL,
                res.par = "Inconclusive"))
  }
  #plot.final <- plot_homo_prob_tetra(pr.hom = pr.hom.4, map = map, ind = ind, parent = parent, title.plot = title.plot)
  problematic.mrk<-setdiff(colnames(pr.hom), colnames(pr.hom.4))
  ## Reestimating missing segments using HMM
  id<-which(dat.mappoly$geno$mrk%in%problematic.mrk)
  if(100*(ncol(pr.hom)-length(id))/ncol(pr.hom) < perc.info){
    pr.hom[]<-0
    return(list(p = "Inconclusive",
         plot.1 = plot.1, 
         plot.final = plot_homo_prob_tetra(pr.hom = pr.hom, map = map, ind = individual, parent, title.plot = title.plot), 
         n.breaks=segment.lengths[[2]], 
         res.par = "Inconclusive"))
  }
  if(length(id) > 0)
    dat.mappoly$geno[id,3:7][]<-1/5
  ## HMM estimation
  w<-calc_genoprob_dist_one_ind(input.map = map.mappoly, dat.dist = dat.mappoly)
  if(parent == "P2"){
    pr.hom.5 <- homo_prob_tetra(ind = individual, genoprob = w, parent = "P2")
  } else {pr.hom.5 <- homo_prob_tetra(ind = individual, genoprob = w, parent = "P1")}
  p<-detect_co_tetra(pr.hom = pr.hom.5, map = map, dist.thresh = dist.thresh)
  plot.final <- plot_homo_prob_tetra(pr.hom = pr.hom.5, map = map, ind = individual, parent = parent, title.plot = title.plot)
  list(p = p$result,
       plot.1 = plot.1, 
       plot.final = plot.final, 
       n.breaks=segment.lengths[[2]], 
       res.par = p$result.par,
       pr.hom.final = pr.hom.5)
}
plot_recombination_points<-function(pr.hom, map, individual, parent, hom.prob.thresh = 0.8, 
                                     seg.length.thresh = 10, perc.info = 20,thresh.nbreak1 = Inf, 
                                    dist.thresh = 2, map.mappoly, dat.mappoly, title.plot = NULL)
{
  p <- count_cross_tetra(pr.hom = pr.hom, 
                    map = map, 
                    individual = individual, 
                    parent = parent, 
                    hom.prob.thresh = hom.prob.thresh, 
                    seg.length.thresh = seg.length.thresh, 
                    dist.thresh = dist.thresh,
                    perc.info = perc.info,
                    map.mappoly = map.mappoly, 
                    dat.mappoly = dat.mappoly,
                   title.plot = title.plot)
  if(as.character(p[[1]][1])=="Inconclusive" | any(p$n.breaks > thresh.nbreak1)){
    return(list(plot=p[[3]], summary = p[[5]], pr.hom.final = p$pr.hom.final, 
                meiotic.configuration = meiotic_configuration(p[[5]], parent = parent), 
                meiotic.graph = meiotic_graph_tetra(a = p[[1]], parent = parent)))
  } else if(as.character(p[[1]][1])=="Zero" | any(p$n.breaks > thresh.nbreak1)){
    return(list(plot=p[[3]], summary = p[[5]], pr.hom.final = p$pr.hom.final,  
                meiotic.configuration = meiotic_configuration(p[[5]], parent = parent), 
                meiotic.graph = meiotic_graph_tetra(a = p[[1]], parent = parent)))
  } else
    {
    p[[1]]<-p[[1]][nchar(as.character(p[[1]]$hom.pair))==2,]
    pb <- ggplot_build(p[[3]])
    pg <- ggplot_gtable(pb)
    for(i in unique(p[[1]]$hom.pair))
    {
      a<-subset(p[[1]], hom.pair%in%i)
      b<-pg$layout[2:7,]
      ystart <- 1/(2*(abs(a$from[1] - a$to[1]) + 1))
      yend <- 1 - ystart
      pnames1 <- names(pg$grobs[[a$from[1]+1]]$children)
      pname1 <- pnames1[str_detect(pnames1, "geom_density")]
      p1 <- pg$grobs[[a$from[1]+1]]$children[[pname1]]$children[[1]]
      x1<-p1$x[1:(length(p1$x)/2)]
      pg <- gtable_add_grob(pg, segmentsGrob(x1[a$point.id], 
                                             rep(ystart, nrow(a)), 
                                             x1[a$point.id],
                                             rep(yend, nrow(a)), 
                                             gp = gpar(lty=1, lwd = 2, col = 1), 
                                             arrow = arrow(angle = 10, length = unit(0.15, "inches"), ends = "both", type = "closed")), 
                            t=b[a$to[1],"t"], b = b[a$from[1],"b"], l=b[a$from[1],"l"])
    }
    z1<-p[[5]]
    colnames(z1)<-c("Hom. Pair.", "CO position (cM)")
    return(list(plot=ggplotify::as.ggplot(pg), summary = z1, pr.hom.final = p$pr.hom.final, plot.orig = p$plot.1, 
                meiotic.configuration = meiotic_configuration(a = z1, parent = parent), 
                meiotic.graph = meiotic_graph_tetra(a = z1, parent = parent)))
  }
}
calc_genoprob_dist_one_ind<-function(input.map, dat.dist, verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  rownames(dat.dist$geno)<-dat.dist$geno$mrk
  mrk.names<-dat.dist$mrk.names[input.map$maps[[1]]$seq.num]
  g <- as.double(t(dat.dist$geno[mrk.names, -c(1:2)]))
  m = dat.dist$m
  n.mrk = input.map$info$n.mrk
  n.ind = 1
  p = as.numeric(unlist(input.map$maps[[1]]$seq.ph$P))
  dp = as.numeric(cumsum(c(0, sapply(input.map$maps[[1]]$seq.ph$P, function(x) sum(length(x))))))
  q = as.numeric(unlist(input.map$maps[[1]]$seq.ph$Q))
  dq = as.numeric(cumsum(c(0, sapply(input.map$maps[[1]]$seq.ph$Q, function(x) sum(length(x))))))
  rf = input.map$maps[[1]]$seq.rf
  indnames<-dat.dist$ind.names
  res.temp <-
    .Call(
      "calc_genoprob_prior",
      as.numeric(m),
      as.numeric(n.mrk),
      as.numeric(n.ind),
      as.numeric(p),
      as.numeric(dp),
      as.numeric(q),
      as.numeric(dq),
      as.double(g),
      as.double(rf),
      as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
      as.double(0),
      as.numeric(FALSE),
      PACKAGE = "mappoly"
    )
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"),
                                mrk.names, indnames)
  structure(list(probs = res.temp[[1]], map = create_map_hap(input.map, dat.dist)), class="mappoly.genoprob")
}
create_map_hap<- function(input.map, dat.dist, step = Inf,
         phase.config = "best")
{
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  mrknames <- dat.dist$mrk.names[input.map$maps[[i.lpc]]$seq.num]
  map <- c(0, cumsum(imf_h(input.map$maps[[i.lpc]]$seq.rf)))
  names(map)<-mrknames
  if(is.null(step))
    return(map)
  minloc <- min(map)
  map <- map-minloc
  a <- seq(floor(min(map)), max(map), by = step)
  a <- a[is.na(match(a,map))]
  names(a) <- paste("loc",a,sep = "_")
  return(sort(c(a,map))+minloc)
}
multi_evidence<-function(res){
  if(as.character(res[1])=="Zero") return(FALSE)
  if(as.character(res[1])=="Inconclusive") return(FALSE)
  a<-res[!duplicated(res$`Hom. Pair.`),,drop = FALSE]
  b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
  a$from <- sapply(b, function(x) x[1])
  a$to <- sapply(b, function(x) x[2])
  a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
  return(any(duplicated(a$from)) | any(duplicated(a$to)))  
}
meiotic_configuration<-function(a, parent){
  if(is.null(a[1])) return(NA)
  if(as.character(a[1])=="Zero") return("z")
  if(as.character(a[1])=="Inconclusive") return("i")
  if(parent == "P1"){
    lt<-letters[1:6]    
  } else {lt<-letters[7:12]}
  #a<-a[!duplicated(a$`Hom. Pair.`),,drop = FALSE]
  b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
  a$from <- sapply(b, function(x) x[1])
  a$to <- sapply(b, function(x) x[2])
  a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
  if(nrow(a) == 1) return("b") 
  x<-a[,3:4]
  g <- igraph::graph_from_data_frame(x)
  igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
  #plot(as.undirected(g))
  I<-igraph::components(g)$csize
  if(all(I < 3)){
    return("b")
  } else if(any(I == 3) | any(I == 4)){
    return("t")} else if(any(I > 4)) {
      return("h")} else return("i") 
}
meiotic_graph_tetra<-function(a, parent){
  if(parent == "P1"){
    lt<-letters[1:4]    
  } else {lt<-letters[5:8]}
  if(grepl(pattern = "a", lt[1])){
    pal<-c("#800000", "#9A6324", "#808000", "#e6194B")
    names(pal)<-lt
  } else if(grepl(pattern = "e", lt[1])){
    pal<-c("#3cb44b", "#469990", "#42d4f4", "#4363d8") 
    names(pal)<-lt
  }
  if(is.null(a[1]))
  {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,6))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=rep(NA, 6)))
  }
  if(as.character(a[1])=="Inconclusive")
  {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,6))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=rep(NA, 6)))
  }  
  if(as.character(a[1])=="Zero") {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,6))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=pal))
  } else {
  b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
  a$from <- sapply(b, function(x) x[1])
  a$to <- sapply(b, function(x) x[2])
  a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
  x<-a[,3:4, drop = FALSE]
  remaining<-setdiff(lt, unique(as.character(as.matrix(x))))
  g <- igraph::graph_from_data_frame(d = x, directed = FALSE)
  if(length(remaining) > 0)
    g <- add_vertices(graph = g, nv = length(remaining), name=remaining)
  V(g)$size<-c(rep(40,6))
  E(g)$edge.color <- "gray80"
  E(g)$width <- 5
  V(g)$label.cex = 2
  igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
  return(list(g=g, vertex.color=pal[names(V(g))]))
  }
}


