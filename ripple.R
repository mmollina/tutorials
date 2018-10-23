TH<-3
cur.map<-MAPs.res[[1]]
win <- 5
tl <- 10 
og<-1:cur.map$info$n.mrk
for(i in 1:(cur.map$info$n.mrk-win+1))
{
  o<-perm_pars(i:(i+win-1))
  l<-numeric(nrow(o))
  for(j in 1:nrow(o))
  { 
    cat(".")
    o1<-c((i-tl):i, o[j,], (win+i):(tl+win))
    x<-get_submap(cur.map, o1[o1>0], verbose = FALSE, use.high.precision = FALSE, tol.final = 10e-3)
    l[j]<-x$maps[[1]]$loglike
  }
  Z<-abs(l-max(l))
  if(Z[1] > TH)
  {
    if(i == 1)
      w<-c(o[which(Z==0),], (win+1):cur.map$info$n.mrk)
    else
      w<-c(1:(i-1), o[which(Z==0),], (win+i):cur.map$info$n.mrk)
    cur.map<-get_submap(cur.map, mrk.pos = w, reestimate.rf = FALSE)
    cat("\n",c((i-tl):i, o[which(Z==0),], (win+i):(tl+win)))
  }
  cat("\n ")
}


b[[2]]<-data.frame(marker = tetra.data$mrk.names[rev(MAPS[[2]]$maps[[1]]$seq.num)],
                   position = cumsum(mappoly::imf_h(c(0, rev(MAPS[[2]]$maps[[1]]$seq.rf)))))
for(i in 1:5)
  b[[i]]$position<-cumsum(c(0,imf_h(MAPs.res[[i]]$maps[[1]]$seq.rf)))
plot(0, type = "n", 
     axes = FALSE, 
     ylim = c(0, max(res[c("polymapR.length", "MAPPoly.length"),])), 
     xlim = c(0, 17), 
     xlab = "Maps",
     ylab = "cM")
axis(2)
xid<-matrix(c(0,1,3,4,6,7,9,10,12,13), nrow = 2)
for(j in 1:5)
{
  lines(x = rep(xid[1,j], 2), y=range(a[[j]]$position))
  lines(x = rep(xid[2,j], 2), y=range(b[[j]]$position))
  idM<-match(a[[j]]$marker, b[[j]]$marker)
  idP<-match(b[[j]]$marker, a[[j]]$marker)
  points(x = rep(xid[1,j], nrow(a[[j]])), y=a[[j]]$position, pch = 19, col = scales::alpha(c("#e41a1c", "#377eb8")[is.na(idM)+1], .2) )
  points(x = rep(xid[2,j], nrow(b[[j]])), y=b[[j]]$position, pch = 19, col = scales::alpha(c("#e41a1c", "#377eb8")[is.na(idP)+1], .2))
  for(i in 1:nrow(a[[j]]))
  {
    if(!is.na(idM[i]))
      lines(c(xid[1,j],xid[2,j]), c(a[[j]]$position[i], b[[j]]$position[idM[i]]), lwd = .5)
  }
}
text(x = apply(xid,2,mean), y = rep(105, 5), labels = c("LG1","LG2","LG3","LG4","LG5"))
legend(14, 50, c("private markers", "shared markers"), 
       col = c("#e41a1c", "#377eb8"), pch = 19)

