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

fam.logcpm <- aggregate(t.logcpm2,by=list(phenos$group),mean)
fam.logcpm <- fam.logcpm[,-1]
scaled.fam.log <- scale(fam.logcpm)

lm.bio.rep.txt.ht <- mclapply(1:56,FUN=function(x){
  select <- as.numeric(which(phenos.192$group== x))
  test.logcpm <- t.logcpm2[-select,]
  test.phenos <- phenos.192[-select,]
  rest <- c(1:56)[-x]
  test.apply <- lapply(1:55,function(x){
    r <- rest[x]
    se <- as.numeric(which(test.phenos$group == r))
    tester <- test.logcpm[-se,]
    tester.p <- test.phenos$ht[-se]
    te.2 <- apply(tester,2, function(x){
      lml <- lm(tester.p ~ x)
      tr <- unlist(summary(lml)[4])
      tr[8]})
    te.2})
  test.apply},mc.cores = 56)


all.55.22.qvalues.ht <- mclapply(1:56,function(x){
  vars <- lm.bio.rep.txt.ht[[x]]
  list.qvalues <- lapply(1:55,function(x){
    var2 <- vars[[x]]
    var.q <- qvalue(var2)
    var.qvalues <- var.q$lfdr
    var.qvalues
  })
  l.qvalues <- unlist(list.qvalues)
  l.qvalues
  #c <- table(names(l.qvalues)[which(l.qvalues <= .015)])                       
  #names(c)[which(c == 55)]
},mc.cores=56) 
hist(table(names(all.55.22.qvalues.ht[[56]])[which(all.55.22.qvalues.ht[[56]] <= .2)]))
se <- seq(.01,.3,.001)
for(each in 1:length(se)){
pred <- vector() ; true <- vector(); 
all.vars <- vector(); l.vars <- vector()
for(i in 1:56){
  c <- table(names(all.55.22.qvalues.ht[[i]])[which(all.55.22.qvalues.ht[[i]] <= .15)])                       
  vars <- unique(c(names(c)[which(c >= 54)],top.txpt.de.fdr.01.ht))
  l.vars <- c(l.vars,length(vars))
 all.vars <- c(all.vars,vars)
  cor.matrix <- cor(t(scaled.fam.log[,match(vars,colnames(scaled.fam.log))])); colnames(cor.matrix) <- 1:56; rownames(cor.matrix) <- 1:56
  #cor.2 <- cor(t(scaled.fam.log[,match(txts.01,colnames(scaled.fam.log))])); rownames(cor.2) <- 1:56; colnames(cor.2) <- 1:56
  test <-  i
  tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.matrix),
                  H2vec = c(1),pheno = phenos.fam,phenoname = "ht")
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])
  #print(cor(tmp[,2],tmp[,3]))
}
plot(pred,true); print(cor(pred,true))
#print(summary(lm(pred~true))[8])
 print(se[[each]])
}
abline(h=mean(true),col="red")
abline(v=mean(pred),col="red")
# best result .15 lfdr weight 1; .26
try <- names(table(all.vars))
length(try)

load("~/Desktop/top.txpt.defdr.01.ht.RDA")
length(which(top.txpt.de.fdr.01.ht %in% try))
try2 <- unique(c(try,top.txpt.de.fdr.01.ht))

pdf(file = "~/Desktop/txt.lmoo.ht.pdf",width = 4,height = 4)
plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
     ylim = c(90,145), xlim=c(95,140))
title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
lml <- summary(lm(pred~true))
legend(x='topleft', legend=paste("R^2 =",round(lml$r.squared,2),sep=""),bty = "n",cex = .75)
dev.off()

list.top.rr.txpts.ht <- vector("list")
for(i in 1:56){
  c <- table(names(all.55.22.qvalues.ht[[i]])[which(all.55.22.qvalues.ht[[i]] <= .15)])                       
  vars <- names(c)[which(c >= 54)]
  list.top.rr.txpts.ht[[i]] <- vars}
save(list.top.rr.txpts.ht,file="~/Desktop/prediction.paper/list.top.rr.txpts.ht.RDA")
