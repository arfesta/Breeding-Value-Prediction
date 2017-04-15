#Load counts from asreml output and remove unneccsary objects####
load("/media/fairfax/seagate/Adam.RNA/analyses/prediction.paper/192.normalization.RDA.RData")
rm(ainv,all.data,ew.countdata,expt.data,full.countdata,lgep.countdata,log.counts.samp,pedigree,practice,tip.countdata,tip.countdata1,
   Bad.samples,i,log.lik.lane0.5,ord,tmp,yo)

# Load libraries####
library(edgeR); library(OmicKriging); library(parallel);library(qvalue)


phenos <- read.table(file="/media/fairfax/seagate/Adam.RNA/analyses/combined/phenos/all.phenos.192.csv", header=T, sep= ",", row.names=1)
phenos <- phenos[order(phenos$animal, decreasing=F),]
phenos$animal[1:48] <- phenos$animal[1:48]-900; phenos$animal[49:192] <- phenos$animal[49:192]-10000
rownames(phenos) <- phenos$animal ;phenos <- phenos[order(phenos$vol, decreasing=F),]
reads.sample.est <- 2^(tip.full.sample.est)
c <- match(rownames(phenos),rownames(reads.sample.est)); reads.sample.est <- reads.sample.est[c,] #set to be same order as phenos object
identical(rownames(phenos),rownames(reads.sample.est))
#TRUE
phenos$group[which(phenos$group==25)] <- 24
phenos$group[which(phenos$group %in% c(26:57))] <- phenos$group[which(phenos$group %in% c(26:57))] -1
phenos.fam <- aggregate(phenos[,4:7],by=list(phenos$group),mean)
grp <- vector("list")
for(i in 1:56){ grp[[i]] <- as.numeric(which(phenos$group == i))}

phenos.192 <- phenos; rownames(phenos.192) <- 1:192

#normalize with EdgeR####
groups = phenos$group # groups dictate the family biological reps belong to
table(groups)

y <- DGEList(counts=t(reads.sample.est),group=groups)
y2 <- calcNormFactors(y); 
logcpm2 <- cpm(y2, log=TRUE, prior.count=5, lib.size=y2$samples[,2])
identical(colnames(logcpm2),rownames(phenos)) 
#TRUE
rm(reads.sample.est,tip.full.sample.est,c,groups,y,y2)


t.logcpm2 <- t(logcpm2)
fam.logcpm2 <- aggregate(t.logcpm2,by=list(phenos$group),mean)
fam.logcpm2 <- fam.logcpm2[,-1]
scaled.fam.logcpm2 <- scale(fam.logcpm2)


cv.folds <- vector("list")
cv.folds[[1]] <- c(1,11,26,35,36,46,54)
cv.folds[[2]] <- c(2,9,12,20,27,37,47)
cv.folds[[3]] <- c(3,13,14,24,28,38,48)
cv.folds[[4]] <- c(4,15,29,39,43,45,49)
cv.folds[[5]] <- c(5,8,16,30,40,50,56)
cv.folds[[6]] <- c(6,17,22,31,33,41,51)
cv.folds[[7]] <- c(7,18,23,32,42,52,55)
cv.folds[[8]] <- c(10,19,21,25,34,44,53)
sort(unlist(cv.folds))

cor.matrix <- cor(t(scaled.fam.logcpm2)); colnames(cor.matrix) <- 1:56; rownames(cor.matrix) <- 1:56

pred <- vector() ; true <- vector()
for(i in 1:8){
    test <-  cv.folds[[i]]
    #test <- i
    tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.matrix),
                    H2vec = c(1),pheno = phenos.fam,phenoname = "ht")
    pred <- c(pred,tmp[,2])
    true <- c(true,tmp[,3])}
  plot(pred,true)
  cor(pred,true)^2
  lml <- summary(lm(pred~true))
  abline(h=mean(true),col="red")
  
  pdf(file = "~/Desktop/txt.all.vol.pdf",width = 4,height = 4)
  plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
       ylim = c(90,145), xlim=c(100,135))
  title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
  lml <- summary(lm(pred~true))
  legend(x='topleft', legend=paste("R^2 =",round(lml$r.squared,2),sep=""),bty = "n",cex = .75)
  dev.off()
  
  pdf(file = "~/Desktop/txt.all.ht.pdf",width = 4,height = 4)
  plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
       ylim = c(95,120), xlim=c(100,115))
  title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
  lml <- summary(lm(pred~true))
  legend(x='topleft', legend=paste("R^2 =",round(lml$r.squared,2),sep=""),bty = "n",cex = .75)
  dev.off()
  
  abline(v=mean(pred),col="red")
  
