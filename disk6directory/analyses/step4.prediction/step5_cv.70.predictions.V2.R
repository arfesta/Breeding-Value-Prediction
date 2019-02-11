library(glmnet); library(doParallel); library(ggplot2)
library(ggpubr); library(OmicKriging); library(caret)
load("./cv.70.data.RData")
names(cv.70.data)

cv.group.predictions.1.70 <- mclapply(1:70, function(rr){
  test.set <- cv.70.data$test.index[rr,]
  training.set <-  c(1:56)[-test.set]
  test.set <- c(1:56)[-training.set]
  
  # First transcripts ridge/lasso
  pval.train <- cv.70.data$all.counts[training.set,]
  pval.test <- cv.70.data$all.counts
  
  #ob.wts <-  1 - (apply(cor(t(pval.train)),2,mean))
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  
#ob.wts <-  1 - (apply(cor(t(pval.train[,gl.mod1])),2,mean))
  #lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,alpha=1,
  #                      nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  
  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,alpha=1,
                        nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  
  txpt.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  #RMSE(txpt.lasso.pred,cv.70.data$fam.phenos$Volume[test.set])
  
  
  #txpt.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
 
  # Next snps ridge/lasso
  pval.train <- as.matrix(cv.70.data$all.snps[training.set,])
  pval.test <- as.matrix(cv.70.data$all.snps)
  #ob.wts <-  1 - (apply(cor(t(pval.train)),2,mean))
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  ob.wts <- 1 - (apply(cor(t(pval.train[,gl.mod1])),2,mean))

  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],
                        parallel = F,alpha=1,nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  snp.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  #cor(snp.lasso.pred,cv.70.data$fam.phenos$Volume[test.set])
  #snp.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
  
 
  both.lasso.pred <- (txpt.lasso.pred + snp.lasso.pred)/2

 cor(both.lasso.pred,cv.70.data$fam.phenos$Volume[test.set])
  
  },mc.cores=10,mc.allow.recursive = F)
hist(unlist(cv.group.predictions.1.70),main="Histogram of 70-fold correlations (r) using txpts + snps",xlab="correlation (r)")
hist(unlist(cv.group.predictions.61.70),main="Histogram of 70-fold correlations (r) using txpts + snps",xlab="correlation (r)")
boxplot(ylim=c(-.3,1),unlist(cv.group.predictions.61.70))
boxplot(ylim=c(-.3,1),unlist(cv.group.predictions.1.70))
summary(unlist(cv.group.predictions.61.70))
#859 seconds

save.image("./step5_cv.70.predictions.v2.RData")

cv.70.group.preds <- c(cv.group.predictions.1.20,cv.group.predictions.21.40,cv.group.predictions.41.60,cv.group.predictions.61.70)
summary(unlist(lapply(cv.70.group.preds,function(x) length(x))))

save(cv.70.group.preds,file="./shared/step5_cv.70.group.preds.v2.RData")

cv.group.predictions.3.70 <- mclapply(1:56, function(rr){
  test.set <- rr
  training.set <-  c(1:56)[-test.set]
  test.set <- c(1:56)[-training.set]
  
  # First transcripts ridge/lasso
  pval.train <- cv.70.data$all.counts[training.set,]
  pval.test <- cv.70.data$all.counts
  
  ob.wts <-  1 - (apply(cor(t(pval.train),method="spearman"),2,mean))
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F,weights = ob.wts,intercept = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  
  ob.wts <-  1 - (apply(cor(t(pval.train[,gl.mod1]),method="spearman"),2,mean))
  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,alpha=1,
                        nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae",weights=ob.wts)
  
  #lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,alpha=1,
  #                      nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  
  txpt.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  RMSE(txpt.lasso.pred,cv.70.data$fam.phenos$Volume[test.set])
  
  
  #txpt.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
  
  # Next snps ridge/lasso
  pval.train <- as.matrix(cv.70.data$all.snps[training.set,])
  pval.test <- as.matrix(cv.70.data$all.snps)
  ob.wts <-  1 - (apply(cor(t(pval.train)),2,mean))
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F,weights=ob.wts)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  ob.wts <- 1 - (apply(cor(t(pval.train[,gl.mod1])),2,mean))
  
  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],
                        parallel = F,alpha=1,nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae",weights=ob.wts)
  snp.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  #cor(snp.lasso.pred,cv.70.data$fam.phenos$Volume[test.set])
  #snp.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
  
  
  both.lasso.pred <- (txpt.lasso.pred + snp.lasso.pred)/2
  
  both.lasso.pred
  
},mc.cores=10,mc.allow.recursive = F)
plot(unlist(cv.group.predictions.3.70),cv.70.data$fam.phenos$Volume)
cor(unlist(cv.group.predictions.3.70),cv.70.data$fam.phenos$Volume)
plot(unlist(cv.group.predictions.3.70) - cv.70.data$fam.phenos$Volume,cv.70.data$fam.phenos$Volume)
