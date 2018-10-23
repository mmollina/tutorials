polymapR_to_mappoly<-function(polymapR.data)
{
  segregating.data <- convert_marker_dosages(dosage_matrix = polymapR.data)
  mappoly.data<-mappoly::hexafake
  mappoly.data$m<-4
  mappoly.data$n.ind<-ncol(segregating.data)-2
  mappoly.data$n.mrk<-nrow(segregating.data)
  mappoly.data$ind.names<-colnames(segregating.data)[-c(1:2)]
  mappoly.data$mrk.names<-rownames(segregating.data)
  mappoly.data$dosage.p<-segregating.data[,1]
  mappoly.data$dosage.q<-segregating.data[,2]
  mappoly.data$sequence<-NA
  mappoly.data$sequence.pos<-NA
  mappoly.data$geno.dose<-segregating.data[,-c(1:2)]  
  mappoly.data$geno.dose[is.na(mappoly.data$geno.dose)]<-5
  mappoly.data$nphen<-0
  mappoly.data$phen<-NULL
  for(i in 1:nrow(mappoly.data$geno.dose))
  {
    mappoly.data$geno.dose[i,mappoly.data$geno.dose[i,]%in%names(table(mappoly.data$geno.dose[i,])[segreg_poly(4, mappoly.data$dosage.p[i], mappoly.data$dosage.q[i])==0])]<-5
  }
  
  
  #######################################
  ## Probability distribution of the genotypes
  #######################################
  myfunc<-function(x, m)
  {
    if(x[1] > m)
      return(mappoly::segreg_poly(m = m, dP = x[2], dQ = x[3]))
    else
    {
      y<-rep(0, m+1)
      y[x[1]+1]<-1
      return(y)
    }
  }
  x<-as.data.frame(as.table(as.matrix(mappoly.data$geno.dose)))
  colnames(x)<-c("mrk", "ind", "dose")
  x$dose.p<-rep(mappoly.data$dosage.p, mappoly.data$n.ind)
  x$dose.q<-rep(mappoly.data$dosage.q, mappoly.data$n.ind)
  y<-t(apply(x[,-c(1:2)], 1, myfunc, m = mappoly.data$m))
  colnames(y)<-0:mappoly.data$m
  z<-cbind(x[,1:2], y)
  mappoly.data$geno<-z
  return(mappoly.data)
}
