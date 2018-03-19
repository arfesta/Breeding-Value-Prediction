##Filter Methods####
library(parallel); library(genefilter)
load("/media/disk6/ARF/RNASEQ/analyses/step2.normalization/dplyr.tech.bio.norm/192.expt.data.Rdata")
load("/media/disk6/ARF/RNASEQ/analyses/step2.normalization/dplyr.tech.bio.norm/sommer.counts.RData")
load("/media/disk6/ARF/RNASEQ/analyses/step2.normalization/asreml_normalization/asreml.counts.RData")
test1 <- scale(asreml.counts,center = T,scale = T)
test2 <- scale(sommer.counts,center=T,scale = T)
all.cnts <- sapply(1:ncol(asreml.counts),function(each.txpt){cor(test1[,each.txpt],test2[,each.txpt])})
boxplot(all.cnts)
these.odd <- which(all.cnts < .95)
summary(apply(asreml.counts[,these.odd],2,sum))
##Create Family mean counts and phenos####
fam.phenos <- data.frame(unique(expt.dat.192[,c("vol","ID")]))
rownames(fam.phenos) <- 1:57

asreml.fam.counts <- lapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% fam.phenos$ID[each.fam])
  if(length(these.rows) == 1){
    2^asreml.counts[these.rows,]
  } else {
    apply(2^asreml.counts[these.rows,],2,mean)
  }
})
asreml.fam.counts <- do.call(rbind,asreml.fam.counts)
rownames(asreml.fam.counts) <- fam.phenos$ID
colnames(asreml.fam.counts) <- colnames(asreml.counts)


sommer.fam.counts <- mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% fam.phenos$ID[each.fam])
  if(length(these.rows) == 1){
    2^sommer.counts[these.rows,]
  } else {
    apply(2^sommer.counts[these.rows,],2,mean)
  }
},mc.cores=20)
sommer.fam.counts <- do.call(rbind,sommer.fam.counts)
rownames(sommer.fam.counts) <- fam.phenos$ID
colnames(sommer.fam.counts) <- colnames(sommer.counts)

a
##Generate general summary stats on txpts####
fam.stats <- vector("list")
fam.stats$txpt.max <- apply(log2(asreml.fam.counts),2,max)
rm.all.negs <- which(fam.stats$txpt.max <= 0)
as.filt.1a <- asreml.fam.counts[,-rm.all.negs]
dim(as.filt.1a)
fam.stats$txpt.sd <- apply(as.filt.1a,2,sd)
rm.low.sd <-(which(round(fam.stats$txpt.sd,digits=2) < .5))
plot(as.filt.1a[,rm.low.sd[8]])
as.filt.2a <- as.filt.1a[,-rm.low.sd]
dim(as.filt.2a)

mypc <- prcomp(log2(as.filt.2a))$x
length(which(abs(mypc[,1]) > 3))
plot(mypc[,1:2])
as.filt.3a <- as.filt.2a[-54,]
as.pheno.3a <- fam.phenos[-54,]; rownames(as.pheno.3a) <- 1:56
mypc <- (prcomp(log2((as.filt.3a)))$x)
plot(mypc[,1:2])

fam.stats$txpt.sd <- apply(as.filt.3a,2,sd)
summary(fam.stats$txpt.sd)
which.max(fam.stats$txpt.sd)
rm.low.sd <-(which(round(fam.stats$txpt.sd,digits=2) < .5))
as.filt.4a <- as.filt.3a[,-rm.low.sd]

dim(as.filt.4a)
#56 66717


###Sepcific filters#####
rm.fam <- which(expt.dat.192$ID == c("N07002"))
bio.filt.dat <- asreml.counts[-rm.fam,which(colnames(asreml.counts) %in% colnames(asreml.fam.counts))]
bio.filt.pheno <- expt.dat.192[-rm.fam,]
#Anova filter pvals of txpt on vol####
fam.txpt.vol.pvals <- mclapply(1:56,function(test.fam){
  rm.rows <- which(bio.filt.pheno$ID %in% as.pheno.3a$ID[test.fam])
  tmp.dat <- bio.filt.dat[-rm.rows,]
  p.vals <- apply(tmp.dat,2,function(each.txpt){
    anova(lm(bio.filt.pheno$vol[-rm.rows] ~ each.txpt))$`Pr(>F)`[1]
  })
},mc.cores=22)

fam.txpt.vol.qvals <- lapply(fam.txpt.vol.pvals,function(each.list){p.adjust(each.list,method="fdr")})
#Anova filter with ID (how well does a model with only family ID explain each txpt?)######
fam.txpt.id.fvals = rowFtests(t(bio.filt.exp2),  as.factor(expt.dat.192$ID))
fam.txpt.id.qvals <- p.adjust(fam.txpt.id.fvals[,2],method="fdr")
length(which(p.adjust(r4[,2],method="bonferroni") < .01))

### # # # # # ####
#NonSpecific Filters####
#CV filter####
fam.txpt.cv.1 <- which(genefilter(t(fam.filt.exp),cv(a=1,b=Inf,na.rm=T)) == T)

##############
#save.image(file = "./step3.1.RData",compress=T)
load("~/step3.1.RData")
###PREDICT####
test.set <- (which(fam.txpt.id.qvals < .000000000000000000000000001))

###GLMNET#####
library(glmnet)

#use cv.filt
tmp.dat <- scale(log2(as.filt.4a[,select]),center=F)

gl.pred <- unlist(lapply(1:56,function(each.fam){
  #tmp.dat <- scale(log2(fam.filt.exp2[,unique(c(which(fam.txpt.vol.qvals[[each.fam]] < .1),fam.txpt.cv.1))]),center=F)
  #tmp.dat <- scale(log2(fam.filt.exp2[,(c(which(fam.txpt.vol.qvals[[each.fam]] < .1)))]),center=F)
  
  gl.mod <- glmnet(x=(tmp.dat[-each.fam,]),y=(as.pheno.3a$vol[-each.fam]),alpha=.99,standardize = F)
  
  the.p <- predict(gl.mod,(tmp.dat),type="response")
  
  the.p[each.fam,ncol(the.p)] - (mean(as.pheno.3a$vol[-each.fam]))
}))

plot(gl.pred,as.pheno.3a$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(gl.pred,as.pheno.3a$vol)^2

###OMICKRIGING####
library(OmicKriging)
cor.dat <- cor(t(scale(log2(as.filt.4a[,select]),center=F)))
#us
ok.pred <- unlist(mclapply(1:56,function(each.fam){
  #cor.dat <- cor(t(scale(log2(fam.filt.exp2[,fam.txpt.vol.qvals[[each.fam]] < .1]),center  =F,scale = T)))
  colnames(cor.dat) <- 1:56; rownames(cor.dat) <- 1:56
  #tmp.dat <- log2(fam.filt.exp2[,fam.txpt.vol.qvals[[each.fam]] < .1])
  okriging(idtest = each.fam,idtrain = c(1:56)[-each.fam],corlist = list(cor.dat),H2vec = c(.99),pheno = as.pheno.3a,phenoname = "vol")[,2]
},mc.cores=20))

plot(ok.pred,as.pheno.3a$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(ok.pred,as.pheno.3a$vol)^2


###randomForest#####
library(randomForest)

#use cv.filt
tmp.dat <- ((bio.filt.exp2[,fam.txpt.cv.1]))

gl.pred <- unlist(mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% fam.phenos$ID[each.fam])
  rf.mod <- randomForest(x = tmp.dat[-these.rows,],y = expt.dat.192$vol[-these.rows],do.trace=F)
  the.imp <- rf.mod$importance[which(rf.mod$importance >=1),1]
  #varImpPlot(rf.mod)
  #plot(mean(predict(rf.mod)[these.rows]),fam.phenos$vol[1])
  #us
  #tmp.dat <- scale(log2(fam.filt.exp2[,test.set]),center=F)
  #both <- unique(c(test.set,fam.txpt.cv.1))
  
  
  tmp.dat <- (log2(fam.filt.exp2[,names(the.imp)]))
  #tmp.dat <- scale(log2(fam.filt.exp2[,unique(c(which(fam.txpt.vol.qvals[[each.fam]] < .1),fam.txpt.cv.1))]),center=F)
  #tmp.dat <- scale(log2(fam.filt.exp2[,(c(which(fam.txpt.vol.qvals[[each.fam]] < .1)))]),center=F)
  
  gl.mod <- glmnet(x=(tmp.dat[-each.fam,]),y=(fam.phenos$vol[-each.fam]),alpha=1,standardize = T)
  
  the.p <- predict(gl.mod,(tmp.dat),type="response")
  
  the.p[each.fam,ncol(the.p)] - (mean(fam.phenos$vol[-each.fam]))
},mc.cores=30))

plot(gl.pred,fam.phenos$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(gl.pred,fam.phenos$vol)^2

library(xgboost)
the.pred <- c()
for(each.fam in 1:57){
  tmp.dat <- xgb.DMatrix((log2(fam.filt.exp2[-each.fam,fam.txpt.vol.qvals[[each.fam]] < .1])), label = fam.phenos$vol[-each.fam])
  #tmp.dat <- cbind("vol"=fam.phenos$vol,log2(fam.filt.exp2[,names(the.imp)]))
  xg.mod <- xgboost(data = tmp.dat,nrounds=50,eta=1,booster="gblinear",verbose = 0)
  tmp.dat2 <- xgb.DMatrix((log2(fam.filt.exp2[,fam.txpt.vol.qvals[[each.fam]] < .1])))
  
  the.pred <- c(the.pred,predict(xg.mod,tmp.dat2)[each.fam])
}
plot(the.pred,fam.phenos$vol)
cor(the.pred,fam.phenos$vol)
