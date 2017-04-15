#Load libraries####
library(OmicKriging)

#Set WD####
setwd("/media/fairfax/sg3/Adam.RNA.SNPs/192.relatedness2")

#First load with SNPeff####
df <- read.table("with.snp.eff/with.missing.snps/out.relatedness2", header=T)
k49 = matrix( 
  df$RELATEDNESS_PHI ,
  nrow=192,
  ncol=192)
rownames(k49) <- unique(df$INDV2)
colnames(k49) <- unique(df$INDV1)
rm(df)

#k22####
df <- read.table("with.snp.eff/with.no.missing.snps/out.relatedness2", header=T)
k22 = matrix( 
  df$RELATEDNESS_PHI ,
  nrow=192,
  ncol=192)
rownames(k22) <- unique(df$INDV2)
colnames(k22) <- unique(df$INDV1)
rm(df)

#summary(c(k22[lower.tri(k22,diag = F)]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05109 0.16450 0.18910 0.18360 0.20280 0.35760 


#k200qual####
df <- read.table("without.snp.eff/with.no.missing.snps/qual.20.filter/out.relatedness2", header=T)
k200.qual = matrix( 
  df$RELATEDNESS_PHI ,
  nrow=192,
  ncol=192)
rownames(k200.qual) <- unique(df$INDV2)
colnames(k200.qual) <- unique(df$INDV1)
rm(df)
#summary(c(k200.qual[lower.tri(k200.qual,diag = F)]))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04549  0.09764  0.13230  0.12420  0.14900  0.32650 

cor(c(k200.qual),c(k22)) #.987

#Both with no missing data are the exact same..creating SNP relationships don't use missing data####
phenos <- read.table(file="/media/fairfax/seagate/Adam.RNA/analyses/combined/phenos/all.phenos.192.csv", header=T, sep= ",", row.names=1)
phenos192 <- phenos[order(phenos$animal, decreasing=F),]
phenos192$animal[1:48] <- phenos192$animal[1:48]-900; phenos192$animal[49:192] <- phenos192$animal[49:192]-10000
phenos192$group[which(phenos192$group==25)] <- 24
phenos192$group[which(phenos192$group %in% c(26:57))] <- phenos192$group[which(phenos192$group %in% c(26:57))] -1

##make groups for biological rep testing####
grps <- vector("list")
for(i in 1:56){grps[[i]] <- which(rownames(phenos192) %in% rownames(phenos192)[which(phenos192$group ==i)])}

rel.mat <- k200.qual; rownames(rel.mat) <- 1:192 ; colnames(rel.mat) <- 1:192
rownames(phenos192) <- 1:192


bio.rep <- 1:192
se <- seq(.01,1,.01)
for(each in 1:length(se)){
pred <- vector(); true <- vector()
for(i in 1:56){
  test <- grps[[i]]
  tmp <- okriging(idtest=test,idtrain=bio.rep[-test],corlist = list(rel.mat),
                  H2vec = c(se[[each]]),pheno = phenos192,phenoname = "ht")
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}
plot(pred,true); cor(pred,true) 
pred.fam <- aggregate(pred,by=list(phenos192$group),mean)
true.fam <- aggregate(true,by=list(phenos192$group),mean)
plot(pred.fam[,2],true.fam[,2]); print(cor(pred.fam[,2],true.fam[,2]))}
###k22 or k49
#vol: weight: 1 cor: .19
#ht: weight: 1 cor: .24
###k200
#vol: weight: 1 cor: .17
#ht: weight: 1 cor: .274
###k200qual
#vol: weight: 1 cor: .18
#ht: weight: 1 cor: .27


####this works#####
##RRBLUP DOESN'T WORK #####
#in order to create correlation matrix from 0,1,2 counts or mean of those there can't be missing data
#because of this we can only use with.snp.eff 22k or without.snp.eff.qual20 filter
phenos.fam <- aggregate(phenos192[,4:6],by=list(phenos192$group),mean)

df <- read.table("without.snp.eff/with.no.missing.snps/qual.20.filter/out.012", header=F,stringsAsFactors = F,row.names=1,na.strings = "-1")
k200.qual.counts <- df
snp.fam <- aggregate(df,by=list(phenos192$group),mean,na.rm=TRUE)
k200.q20.snp.fam.df <- snp.fam[,-1]

df <- read.table("with.snp.eff/with.no.missing.snps/out.012", header=F,stringsAsFactors = F,row.names=1,na.strings = "-1")
k22.counts <- df
rm(df)
k22.counts <- k22.counts-1


snp.fam <- aggregate(k22.counts,by=list(phenos192$group),mean,na.rm=TRUE)
k22.snp.fam.df <- snp.fam[,-1]
phenos.fam <- aggregate(phenos192[,4:7],by=list(phenos192$group),mean)
library(rrBLUP)
A <- A.mat(k22.snp.fam.df); rownames(A) <- 1:56; colnames(A) <- 1:56
pred <- vector(); true <- vector()
for(i in 1:56){
k22.pred <- data.frame(y=phenos.fam$vol,gid=1:56)
ls <- i
k22.pred[ls,1] <- NA
c <- kin.blup(data = k22.pred,geno="gid",pheno="y",K = A)
pred <- c(pred,c$pred[ls])
true <- c(true,phenos.fam$vol[ls])
print(i)}
cor(pred,true)
pre.fam <- aggregate(pred,by=list(phenos192$group),mean)
cor(pre.fam[,2],true.fam[,2])
