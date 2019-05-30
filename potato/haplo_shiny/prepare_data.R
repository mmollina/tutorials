z<-readRDS("~/Documents/EiB_May_2019/Demos/mapping/potato/genoprob_list.rds")
source("~/repos/haplo_shiny/detect_meiotic_configuration_functions.R")
for(ch in 1:length(z)){
  cat("chromosome: ", ch, "\n")
  w<-z[[ch]]
  ind<-dimnames(w$probs)[[3]]
  zt<-zb<-array(data = NA, dim = c(dim(w$probs)[3], 4, dim(w$probs)[2]))
  dimnames(zb)<-list(ind, letters[1:4], dimnames(w$probs)[[2]])
  dimnames(zt)<-list(ind, letters[5:8], dimnames(w$probs)[[2]])
  for(i in 1:length(ind))
  {
    zb[i,,]<-homo_prob_tetra(ind = ind[i], genoprob = w, parent = "P1")
    zt[i,,]<-homo_prob_tetra(ind = ind[i], genoprob = w, parent = "P2")
  }
  map<-w$map
  save(zb, zt, map, file = paste0("hom_prob_ch", ch, ".rda"))
}
ind.names<-dimnames(w$probs)[[3]]
save(ind.names, file = "individuals.rda")

final.map <- readRDS("~/Documents/EiB_May_2019/Demos/mapping/potato/final.map.rds")
for(ch in 1:length(z)){
  map.mappoly = final.map[[ch]]
  save(map.mappoly, file = paste0("./dat/map_ch",ch,".rda"))
}

input.data<-readRDS("~/Documents/EiB_May_2019/Demos/mapping/potato/input.data.rds")
for(individual in ind.names){
  cat("ind: ", individual, "\n")
  dat.mappoly<-input.data
  dat.mappoly$geno<-subset(dat.mappoly$geno, ind%in%individual)
  dat.mappoly$n.ind<-length(individual)
  dat.mappoly$ind.names<-individual
  dat.mappoly$geno.dose<-dat.mappoly$geno.dose[,individual, drop = FALSE]
  save(dat.mappoly, file = paste0("./dat/dat_ind_",individual,".rda"))
}
