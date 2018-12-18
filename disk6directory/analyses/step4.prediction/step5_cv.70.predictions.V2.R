library(glmnet); library(doParallel); library(ggplot2)
library(ggpubr); library(OmicKriging); library(caret)
load("./cv.70.data.RData")
names(cv.70.data)

pt <-proc.time()
cv.group.predictions.61.70 <- mclapply(61:70, function(rr){
  test.set <- cv.70.data$test.index[rr,]
  training.set <-  c(1:56)[-test.set]
  test.set <- c(1:56)[-training.set]
  
  # First transcripts ridge/lasso
  pval.train <- cv.70.data$all.counts[training.set,]
  pval.test <- cv.70.data$all.counts
  
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  
  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,alpha=1,
                        nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  
  txpt.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  txpt.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
  ridge.cv <- cv.glmnet(x = pval.train,y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,
                        alpha=0,nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  txpt.ridge.pred <- predict.cv.glmnet(object = ridge.cv,newx = pval.test,s = ridge.cv$lambda.min)[test.set]
  
  # transcripts OK
  x.train <- cor(t(pval.train))
  y.train <- cv.70.data$fam.phenos[training.set,]
  all.weights <- seq(.1,1,.1)
  
  train.rmse <- unlist(lapply(1:length(all.weights),function(each.weight){
    this.w <- all.weights[each.weight]
    ok.preds <- lapply(1:nrow(y.train),function(each.cv){
      ok.pred <- okriging(idtest = as.character(y.train$fam_id)[each.cv],idtrain = as.character(y.train$fam_id)[-each.cv],corlist = list(x.train),H2vec = c(this.w),pheno = y.train,phenoname = "Volume")
    })
    ok.preds <-do.call(rbind,ok.preds)
    RMSE(ok.preds[,2],ok.preds[,3])
  }))
  
  final.weight <- all.weights[which.min(train.rmse)]
  all.d <- cv.70.data$all.counts
  all.p <- (cv.70.data$fam.phenos)
  x.cor <- cor(t(all.d))
  ok.pred <- okriging(idtest = as.character(cv.70.data$fam.phenos$fam_id[test.set]),
                      idtrain =as.character(cv.70.data$fam.phenos$fam_id[-c(test.set)]),
                      corlist = list(x.cor),H2vec = c(final.weight),pheno = all.p ,phenoname = "Volume")
  txpt.ok.pred <- ok.pred[,-c(1)]
  
  # Next snps ridge/lasso
  pval.train <- as.matrix(cv.70.data$all.snps[training.set,])
  pval.test <- as.matrix(cv.70.data$all.snps)
  
  gl.mod <- glmnet(x=pval.train,y=cv.70.data$fam.phenos$Volume[training.set],alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- which(apply(gl.mod1,1,function(x) length(unique(x)))[-1] > 1)
  
  lasso.cv <- cv.glmnet(x = pval.train[,gl.mod1],y = cv.70.data$fam.phenos$Volume[training.set],
                        parallel = F,alpha=1,nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  snp.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,gl.mod1],s = lasso.cv$lambda.min)[test.set]
  snp.lasso.features <- colnames(pval.train[,gl.mod1])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
  ridge.cv <- cv.glmnet(x = pval.train,y = cv.70.data$fam.phenos$Volume[training.set],parallel = F,
                        alpha=0,nfolds = length(training.set),grouped = F,standardize=F,type.measure="mae")
  snp.ridge.pred <- predict.cv.glmnet(object = ridge.cv,newx = pval.test,s = ridge.cv$lambda.min)[test.set]
  
  # snps OK
  x.train <- cor(t(pval.train))
  y.train <- cv.70.data$fam.phenos[training.set,]
  all.weights <- seq(.1,1,.1)
  
  train.rmse <- unlist(lapply(1:length(all.weights),function(each.weight){
    this.w <- all.weights[each.weight]
    ok.preds <- lapply(1:nrow(y.train),function(each.cv){
      ok.pred <- okriging(idtest = as.character(y.train$fam_id)[each.cv],idtrain = as.character(y.train$fam_id)[-each.cv],corlist = list(x.train),H2vec = c(this.w),pheno = y.train,phenoname = "Volume")
    })
    ok.preds <-do.call(rbind,ok.preds)
    RMSE(ok.preds[,2],ok.preds[,3])
  }))
  
  final.weight <- all.weights[which.min(train.rmse)]
  all.d <- cv.70.data$all.snps
  all.p <- (cv.70.data$fam.phenos)
  x.cor <- cor(t(all.d))
  ok.pred <- okriging(idtest = as.character(cv.70.data$fam.phenos$fam_id[test.set]),idtrain =as.character(cv.70.data$fam.phenos$fam_id[-c(test.set)]),corlist = list(x.cor),H2vec = c(final.weight),pheno = all.p ,phenoname = "Volume")
  snp.ok.pred <- ok.pred[,-1]
  
  both.ridge.pred <- (txpt.ridge.pred + snp.ridge.pred)/2
  both.lasso.pred <- (txpt.lasso.pred + snp.lasso.pred)/2
  both.ok.pred <- (txpt.ok.pred + snp.ok.pred)/2
  
  prediction.ridge <- data.frame("txpts"=txpt.ridge.pred,"snps"=snp.ridge.pred,"txpts.snps"=both.ridge.pred)
  colnames(prediction.ridge) <- c("txpts","snps","txpts.snps")
  prediction.lasso <- data.frame("txpts"=txpt.lasso.pred,"snps"=snp.lasso.pred,"txpts.snps"=both.lasso.pred)
  colnames(prediction.lasso) <- c("txpts","snps","txpts.snps")
  prediction.ok <- data.frame("txpts"=txpt.ok.pred[,-2],"snps"=snp.ok.pred[,-2],"txpts.snps"=both.ok.pred[,-2])
  colnames(prediction.ok) <- c("txpts","snps","txpts.snps")
  
  lasso.features <- list("txpts"=txpt.lasso.features,"snps"=snp.lasso.features)
  
  results <- list("ridge"=prediction.ridge,"lasso"=prediction.lasso,"ok"=prediction.ok,"lasso.feats"=lasso.features)
},mc.cores=10,mc.allow.recursive = F)
proc.time()-pt
#859 seconds

save.image("./step5_cv.70.predictions.v2.RData")

cv.70.group.preds <- c(cv.group.predictions.1.20,cv.group.predictions.21.40,cv.group.predictions.41.60,cv.group.predictions.61.70)
summary(unlist(lapply(cv.70.group.preds,function(x) length(x))))

save(cv.70.group.preds,file="./shared/step5_cv.70.group.preds.v2.RData")

