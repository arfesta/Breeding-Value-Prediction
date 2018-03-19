snp.bio.1 <- read.table(file = "/media/disk6/ARF/RNASEQ/snps/86k/out.012",row.names = 1)
snp.pos <- read.table(file = "/media/disk6/ARF/RNASEQ/snps/86k/out.012.pos")
snp.ind <- read.table(file = "/media/disk6/ARF/RNASEQ/snps/86k/out.012.indv")
snp.ind <- gsub(x=gsub("\\..*","", snp.ind[,1] ),pattern = "_",replacement = "")
snp.ind[1:48] <- as.numeric(snp.ind[1:48]) + 100
library(stringr)
regexp <- "[[:digit:]]+"

# process string
snp.ind[49:192] <- as.numeric(str_extract(snp.ind[49:192], regexp)) + 1000

rownames(snp.bio.1) <- as.character(snp.ind)
colnames(snp.bio.1) <- paste(snp.pos$V1,snp.pos$V2,sep="_")
identical(rownames(snp.bio.1)[match(expt.dat.192$animal.id,rownames(snp.bio.1))],expt.dat.192$animal.id)

snp.bio.1 <- snp.bio.1[match(expt.dat.192$animal.id,rownames(snp.bio.1)),]
fam.phenos <- data.frame(unique(expt.dat.192[,c("vol","ID")]))
rownames(fam.phenos) <- 1:57

fam.snp.1 <- mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% fam.phenos$ID[each.fam])
  if(length(these.rows) == 1){
    snp.bio.1[these.rows,]
  } else {
    apply(snp.bio.1[these.rows,],2,mean)
  }
},mc.cores =20)
fam.snp.1 <- do.call(rbind,fam.snp.1)
rownames(fam.snp.1) <- fam.phenos$ID
colnames(fam.snp.1) <- colnames(snp.bio.1)


all.sd <- apply(fam.snp.1,2,sd)
all.vars <- apply(fam.snp.1,2,var)
all.mean <- apply(fam.snp.1,2,mean)
all.dispersion <- all.vars/all.mean
all.cv <- all.sd/all.mean
these.snps <- which(all.dispersion > 1 & all.cv > 1)
filt.snps <- fam.snp.1[,these.snps]
summary(apply(filt.snps,2,var))
boxplot(all.vars)
which.max(all.vars)
#est_Pita_61670413_7

fam.snp.1[,"est_pita_70977076_292"]
library(caret)
rm.txpts <- findCorrelation(x = cor(filt.snps))
filt.snps2 <- filt.snps[,-rm.txpts]
library(OmicKriging)
cor.dat <- cor(t(scale(filt.snps2,center=F,scale=T)))
cor.dat1 <- cor(t(scale(filt.snps2,center=F,scale=T)))
cor.dat2 <- cor(t(scale(log2(fam.filt.exp2[,fam.txpt.cv.1 >1]),center=F)))
these.snps <- which(order(all.cv,decreasing = T) %in% c(1:50))
filt.snps2 <- fam.snp.1[,these.snps]
cor.dat1 <- cor(t(scale(filt.snps2,center=F,scale=T)))
ok.pred <- unlist(mclapply(1:57,function(each.fam){

  colnames(cor.dat1) <- 1:57; rownames(cor.dat1) <- 1:57
  colnames(cor.dat2) <- 1:57; rownames(cor.dat2) <- 1:57
  okriging(idtest = each.fam,idtrain = c(1:57)[-each.fam],corlist = list(cor.dat1),H2vec = c(1),pheno = fam.phenos,phenoname = "vol")[,2]
},mc.cores=20))

plot(ok.pred,fam.phenos$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(ok.pred,fam.phenos$vol)^2

all.sd.bio <- app
these.snps <- which(order(all.cv,decreasing = T) %in% c(1:20000))
filt.snps2 <- as.matrix(fam.snp.1[,these.snps])
rm.snps <- which(apply(filt.snps2,2,var) < .01)
filt.snps2 <- filt.snps2[,-rm.snps]
rm.snps <- findCorrelation(x = cor(filt.snps2),cutoff = .5)
filt.snps2 <- filt.snps2[,-rm.snps]
filt.bio2 <- snp.bio.1[,these.snps]

library(glmnet)
gl.pred <- unlist(mclapply(1:57,function(each.fam){
  #tmp.dat <- scale(log2(fam.filt.exp2[,unique(c(which(fam.txpt.vol.qvals[[each.fam]] < .1),fam.txpt.cv.1))]),center=F)
  #tmp.dat2 <- scale(as.matrix(cbind(filt.snps2,tmp.dat)),center=F,scale=T)
  tmp.dat2 <- scale(as.matrix(filt.snps2),center=F)
  #these.rows <- which(expt.dat.192$ID %in% c(fam.phenos$ID[each.fam]))
  gl.mod <- glmnet(x=(tmp.dat2[-each.fam,]),y=(fam.phenos$vol[-each.fam]),alpha=0,
                   standardize.response = T,standardize = F)
  
  the.p <- predict(gl.mod,(tmp.dat2),type="response")
  
  (the.p[each.fam,ncol(the.p)]) - (mean(fam.phenos$vol[-each.fam]))
},mc.cores=30))
plot(gl.pred,fam.phenos$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(gl.pred,fam.phenos$vol)^2
