require(ClusterCall)
require(polymap2)
setwd("~/repos/potato_analysis_to_publish/cluster_snp_call/")
datxy<-read.table(file = "~/repos/potato_analysis_to_publish/cluster_snp_call/potato_x_y_solcap_public.txt",
                  header=TRUE)
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
write.csv(x = dat.theta, file = "~/repos/potato_analysis_to_publish/cluster_snp_call/solcap_public_theta.csv")

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
write.csv(x = dat.r, file = "~/repos/potato_analysis_to_publish/cluster_snp_call/solcap_public_r.csv")

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


#######################################################################################
## Using boundaries provided by SolaCap Consortium

# bound<-read.csv("SolCAP_Infinium_8300_5_cluster_positions.csv")
# head(bound)
# d0<-which(dat.theta[1,] >= bound[1, 2] & dat.theta[1,] <= bound[1, 3])
# d1<-which(dat.theta[1,] >= bound[1, 4] & dat.theta[1,] <= bound[1, 5])
# d2<-which(dat.theta[1,] >= bound[1, 6] & dat.theta[1,] <= bound[1, 7])
# d3<-which(dat.theta[1,] >= bound[1, 8] & dat.theta[1,] <= bound[1, 9])
# d4<-which(dat.theta[1,] >= bound[1, 10] & dat.theta[1,] <= bound[1, 11])







