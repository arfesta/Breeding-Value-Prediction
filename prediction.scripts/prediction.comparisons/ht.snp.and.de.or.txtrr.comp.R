load("~/Desktop/prediction.paper/list.top.rr.txpts.ht.RDA")
load("~/Desktop/prediction.paper/k22.family.counts.RDA")
load("~/Desktop/prediction.paper/phenos.family.RDA")
load("~/Desktop/prediction.paper/3.2.2/top.txpt.defdr.01.ht.RDA")
load("~/Desktop/prediction.paper/scaled.family.counts.RDA")
load("~/Desktop/prediction.paper/top.roundrobin.snp.ht.RDA")

library(OmicKriging)

#Compare weighing SNPs found from RR LOO analysis vs. top DE transcripts found####
cor.txt <- cor(t(scaled.fam.logcpm2[,top.txpt.de.fdr.01.ht])); colnames(cor.txt) <- 1:56; rownames(cor.txt) <- 1:56
se <- seq(.01,.99,.01)
cor.list3 <- vector()
for(each in 1:length(se)){
  pred <- vector() ; true <- vector()
  for(i in 1:56){
    cor.snp <- cor(t(k22.snp.fam.df[,top.roundrobin.snp.ht[[i]]])); colnames(cor.snp) <- 1:56; rownames(cor.snp) <- 1:56
    test <-  i
    tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.txt,cor.snp),
                    H2vec = c(se[[each]],1-se[[each]]),pheno = phenos.fam,phenoname = "ht")
    pred <- c(pred,tmp[,2])
    true <- c(true,tmp[,3])
    cor(tmp[,2],tmp[,3])}
  plot(pred,true)
  cor.list3 <- c(cor.list3,cor(pred,true)^2)}

plot(se,cor.list3)

#Compare weighing SNPs found from RR LOO analysis vs. LOO transcripts found####
se <- seq(.01,.98,.01)
cor.list4 <- vector()
for(each in 1:length(se)){
  pred <- vector() ; true <- vector()
  for(i in 1:56){
    cor.txt <- cor(t(scaled.fam.logcpm2[,list.top.rr.txpts.ht[[i]]])); colnames(cor.txt) <- 1:56; rownames(cor.txt) <- 1:56
    cor.snp <- cor(t(k22.snp.fam.df[,top.roundrobin.snp.ht[[i]]])); colnames(cor.snp) <- 1:56; rownames(cor.snp) <- 1:56
    test <-  i
    tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.txt,cor.snp),
                    H2vec = c(se[[each]],.99-se[[each]]),pheno = phenos.fam,phenoname = "ht")
    pred <- c(pred,tmp[,2])
    true <- c(true,tmp[,3])
    cor(tmp[,2],tmp[,3])}
  plot(pred,true)
  cor.list4 <- c(cor.list4,cor(pred,true)^2)}

plot(se,cor.list4)
summary(cor.list4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3728  0.4049  0.4434  0.4383  0.4752  0.4883 

pred <- vector() ; true <- vector()
for(i in 1:56){
  #cor.txt <- cor(t(scaled.fam.logcpm2[,list.top.rr.txpts.vol[[i]]])); colnames(cor.txt) <- 1:56; rownames(cor.txt) <- 1:56
  cor.snp <- cor(t(k22.snp.fam.df[,snp.vol.list.vars[[i]]])); colnames(cor.snp) <- 1:56; rownames(cor.snp) <- 1:56
  test <-  i
  tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.snp),
                  H2vec = c(.9),pheno = phenos.fam,phenoname = "vol")
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])
  cor(tmp[,2],tmp[,3])}
plot(pred,true); cor(pred,true)^2
