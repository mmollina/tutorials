require(netgwas)
setwd("~/repos/tutorials/comparing_methods/")
load(file = "~/repos/BT_map/figures/RData/all.maps.rda")
load("~/repos/BT_map/src/identifying_low_quality_individuals/BT_trifida_triloba_298.RData")
x1<-all.maps[[1]][[4]]$maps[[1]]$seq.num
x2<-all.maps[[2]][[4]]$maps[[1]]$seq.num

o<-c(x1[ceiling(seq(1, length(x1), length.out = 500))],
     x2[ceiling(seq(1, length(x2), length.out = 500))])

dat.hexa<-BT.trifida.triloba.298$geno.dose[o,]
dat.hexa[dat.hexa==7]<-NA
dat.hexa<-t(dat.hexa)
dim(dat.hexa)
hexa.map <- netmap(dat.hexa, cross = "outbred")
save(hexa.map, file = "hexa.map.RData")

hexa.map$map

plot(match(as.character(hexa.map$map$markers), BT.trifida.triloba.298$mrk.names[o]))


x<-hexa.map
opt.theta <- x$allres$Theta[[x$opt.index]]
rownames(opt.theta)  <- colnames(opt.theta) 
orderedm <- opt.theta[colnames(x$allres$data)  , colnames(x$allres$data) ]
image(as.matrix(orderedm),  xlab="markers", ylab="markers", main="Ordered markers", cex=1, sub="")







