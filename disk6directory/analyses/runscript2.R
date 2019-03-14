library(glmnet); library(doParallel); library(ggplot2)
library(OmicKriging); library(caret)
library(ggpubr)
load("../RNA_SEQ_Data/shared/lgep.vs.ew.data.RData")
names(lgep.vs.ew.data)

{
# transcripts OK
  pval.train <- lgep.vs.ew.data$train.counts
  pval.test <- lgep.vs.ew.data$test.counts
  x.train <- cor(t(pval.train))
  y.train <- lgep.vs.ew.data$train.phenos
  all.weights <- seq(.1,1,.1)

train.mae <- unlist(lapply(1:length(all.weights),function(each.weight){
  this.w <- all.weights[each.weight]
  ok.preds <- lapply(1:nrow(y.train),function(each.cv){
    ok.pred <- okriging(idtest = y.train$fam_id[each.cv],idtrain = y.train$fam_id[-each.cv],corlist = list(x.train),H2vec = c(this.w),pheno = y.train,phenoname = "Volume")
  })
  ok.preds <-do.call(rbind,ok.preds)
  mean(abs(ok.preds[,2]-ok.preds[,3]))
}))

final.weight <- all.weights[which.min(train.mae)]
all.d <- rbind(pval.train,pval.test)
all.p <- rbind(y.train,lgep.vs.ew.data$test.phenos)
x.cor <- cor(t(all.d))
ok.pred <- okriging(idtest = lgep.vs.ew.data$test.phenos$fam_id,idtrain = lgep.vs.ew.data$train.phenos$fam_id,corlist = list(x.cor),H2vec = c(final.weight),pheno = all.p ,phenoname = "Volume")
txpt.ok.pred <- ok.pred[,-c(1)]

# Next snps ridge/lasso
pval.train <- as.matrix(lgep.vs.ew.data$train.snps)
pval.test <- as.matrix(lgep.vs.ew.data$test.snps)

# snps OK
x.train <- cor(t(pval.train))
y.train <- lgep.vs.ew.data$train.phenos
all.weights <- seq(.1,1,.1)

train.mae <- unlist(lapply(1:length(all.weights),function(each.weight){
  this.w <- all.weights[each.weight]
  ok.preds <- lapply(1:nrow(y.train),function(each.cv){
    ok.pred <- okriging(idtest = y.train$fam_id[each.cv],idtrain = y.train$fam_id[-each.cv],corlist = list(x.train),H2vec = c(this.w),pheno = y.train,phenoname = "Volume")
  })
  ok.preds <-do.call(rbind,ok.preds)
  mean(abs(ok.preds[,2]-ok.preds[,3]))
}))

final.weight <- all.weights[which.min(train.mae)]
all.d <- rbind(pval.train,pval.test)
all.p <- rbind(y.train,lgep.vs.ew.data$test.phenos)
x.cor <- cor(t(all.d))
ok.pred <- okriging(idtest = lgep.vs.ew.data$test.phenos$fam_id,idtrain = lgep.vs.ew.data$train.phenos$fam_id,corlist = list(x.cor),H2vec = c(final.weight),pheno = all.p ,phenoname = "Volume")
snp.ok.pred <- ok.pred[,-1]

both.ok.pred <- (snp.ok.pred + txpt.ok.pred)/2

x.txpt.train <- cor(t(lgep.vs.ew.data$train.counts))
x.snp.train <- cor(t(as.matrix(lgep.vs.ew.data$train.snps)))

# OK
y.train <- lgep.vs.ew.data$train.phenos
all.weights <- expand.grid(seq(.1,1,.1),seq(.1,1,.1))
all.weights <- all.weights[-which(apply(all.weights,1,sum) > 1),]

train.mae <- unlist(lapply(1:nrow(all.weights),function(each.weight){
  this.w <- all.weights[each.weight,1]
  this.w2 <- all.weights[each.weight,2]
  ok.preds <- lapply(1:nrow(y.train),function(each.cv){
    ok.pred <- okriging(idtest = y.train$fam_id[each.cv],idtrain = y.train$fam_id[-each.cv],corlist = list(x.txpt.train,x.snp.train),H2vec = c(this.w,this.w2),pheno = y.train,phenoname = "Volume")
  })
  ok.preds <-do.call(rbind,ok.preds)
  (MAE(ok.preds[,2],ok.preds[,3]))
}))


final.weight <- all.weights[which.min(train.mae),]
x.cor.txpt <- cor(t(rbind(lgep.vs.ew.data$train.counts,lgep.vs.ew.data$test.counts)))
x.cor.snp <- cor(t(rbind(lgep.vs.ew.data$train.snps,lgep.vs.ew.data$test.snps)))

all.p <- rbind(y.train,lgep.vs.ew.data$test.phenos)

ok.pred <- okriging(idtest = lgep.vs.ew.data$test.phenos$fam_id,idtrain = lgep.vs.ew.data$train.phenos$fam_id,corlist = list(x.cor.txpt,x.cor.snp),H2vec = c(final.weight[1,1],final.weight[1,2]),pheno = all.p ,phenoname = "Volume")
txpt.snp.ok.pred <- ok.pred[,-c(1)]
prediction.ok <- data.frame("txpts"=txpt.ok.pred[,-2],"snps"=snp.ok.pred[,-2],"txpt.snp.ok"=txpt.snp.ok.pred[,1],"txpts.snps"=both.ok.pred[,-2])
}
{
# Lasso
  #txpts
  pval.train <- lgep.vs.ew.data$train.counts
  pval.test <- lgep.vs.ew.data$test.counts
  
  gl.mod <- glmnet(x = pval.train,y =lgep.vs.ew.data$train.phenos$Volume,alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef(gl.mod))
  gl.mod1 <- apply(gl.mod1,1,function(x) length(unique(x)))[-1]
  these.txpts <-  which(gl.mod1 > 1)
  
  lasso.cv <- cv.glmnet(x = pval.train[,these.txpts],y = lgep.vs.ew.data$train.phenos$Volume,parallel = F,alpha=1,nfolds = 10,grouped = F,standardize=F,type.measure = "mae")
  
  txpt.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,these.txpts],s = lasso.cv$lambda.min)
  txpt.lasso.features <- colnames(pval.train[,these.txpts])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min) != 0)]
  
 #snps
  pval.train <- as.matrix(lgep.vs.ew.data$train.snps)
  pval.test <- as.matrix(lgep.vs.ew.data$test.snps)
  
  gl.mod <- glmnet(x = pval.train,y =lgep.vs.ew.data$train.phenos$Volume,alpha = 1,standardize = F)
  gl.mod1 <- as.matrix(coef.glmnet(gl.mod))
  gl.mod1 <- apply(gl.mod1,1,function(x) length(unique(x)))[-1]
  these.snps <-  which(gl.mod1 > 1)
  lasso.cv <- cv.glmnet(x = pval.train[,names(these.snps)],y = lgep.vs.ew.data$train.phenos$Volume,parallel = F,alpha=1,nfolds = 10,grouped = F,standardize=T,type.measure = "mae")
  
  snp.lasso.pred <- predict.cv.glmnet(object = lasso.cv,newx = pval.test[,names(these.snps)],s = lasso.cv$lambda.min)
  snp.lasso.features <- colnames(pval.train[,names(these.snps)])[which(coef.cv.glmnet(lasso.cv,s = lasso.cv$lambda.min)[-1] != 0)]
  
  
  x.txpt.train <- cor(t(lgep.vs.ew.data$train.counts[,these.txpts]))
  x.snp.train <- cor(t(as.matrix(lgep.vs.ew.data$train.snps)[,these.snps]))
 
  # OK
  y.train <- lgep.vs.ew.data$train.phenos
  all.weights <- expand.grid(seq(.1,1,.1),seq(.1,1,.1))
  all.weights <- all.weights[-which(apply(all.weights,1,sum) > 1),]
  
  train.mae <- unlist(lapply(1:nrow(all.weights),function(each.weight){
    this.w <- all.weights[each.weight,1]
    this.w2 <- all.weights[each.weight,2]
    ok.preds <- lapply(1:nrow(y.train),function(each.cv){
      ok.pred <- okriging(idtest = y.train$fam_id[each.cv],idtrain = y.train$fam_id[-each.cv],corlist = list(x.txpt.train,x.snp.train),H2vec = c(this.w,this.w2),pheno = y.train,phenoname = "Volume")
    })
    ok.preds <-do.call(rbind,ok.preds)
    (MAE(pred = ok.preds[,2],obs = ok.preds[,3]))
  }))

  
  final.weight <- all.weights[which.min(train.mae),]
  x.cor.txpt <- cor(t(rbind(lgep.vs.ew.data$train.counts[,these.txpts],lgep.vs.ew.data$test.counts[,these.txpts])))
  x.cor.snp <- cor(t(rbind(lgep.vs.ew.data$train.snps[,these.snps],lgep.vs.ew.data$test.snps[,these.snps])))
  
  all.p <- rbind(y.train,lgep.vs.ew.data$test.phenos)

  ok.pred <- okriging(idtest = lgep.vs.ew.data$test.phenos$fam_id,idtrain = lgep.vs.ew.data$train.phenos$fam_id,corlist = list(x.cor.txpt,x.cor.snp),H2vec = c(final.weight[1,1],final.weight[1,2]),pheno = all.p ,phenoname = "Volume")
  txpt.snp.lasso.ok.pred <- ok.pred[,-c(1)]
  
  
  both.lasso.pred <- (txpt.lasso.pred + snp.lasso.pred)/2

  prediction.lasso <- data.frame("txpts"=txpt.lasso.pred,"snps"=snp.lasso.pred,"txpts.snps.lasso.ok"=txpt.snp.lasso.ok.pred[,1],"txpts.snps"=both.lasso.pred)
  #colnames(prediction.lasso) <- c("txpts.lasso","snps.lasso","txpts.snps.lasso.mean","txpts.snps.lasso.ok")
    
  lasso.features <- list("txpts"=txpt.lasso.features,"snps"=snp.lasso.features)
}    
  results <- list("lasso"=prediction.lasso,"ok"=prediction.ok,"lasso.feats"=lasso.features)
res <-  (rbind(
  "OmicKriging"=apply(results$ok,2,function(x) cor(x,lgep.vs.ew.data$test.phenos$Volume)^2),
  "Lasso"=apply(results$lasso,2,function(x) cor(x,lgep.vs.ew.data$test.phenos$Volume)^2)
  ))
colnames(res) <- c("Transcripts","SNPs","OKTS","MTS") 

res2 <-  (rbind(
  "OmicKriging"=apply(results$ok,2,function(x) RMSE(x,lgep.vs.ew.data$test.phenos$Volume)),
  "Lasso"=apply(results$lasso,2,function(x) RMSE(x,lgep.vs.ew.data$test.phenos$Volume))
))
colnames(res2) <- c("Transcripts","SNPs","OKTS","MTS") 

res3 <- rbind("OmicKriging"=paste0(round(res[1,],digits=5)," (",round(res2[1,],digits = 3),")"),
      "Lasso"=paste0(round(res[2,],digits=5)," (",round(res2[2,],digits = 3),")")
)
colnames(res3) <- c("Transcripts","SNPs","OKTS","MTS") 
  res3<- t(res3)
  
  
  # weights used
  #.9t .1e transcripts
  # 1s 0e snps
  # .8t .2s txpts + snps
  # .8t .2s
  #Prediction accuracies of Batch 2 families (n=11) with OmicKriging & Lasso utilizing 
  ##Transcripts, SNPs, both Transcripts and SNPs in sepearate correlation matrices (OKTS), or the mean of the prediction from Transcripts and SNPs (MTS).
  
  # variables identified for prediction
  lapply(results$lasso.feats,length)
  #49 txpts
  #39 snps
  lasso.txpts <- as.character(results$lasso.feats$txpts)
  rm(load.counts_tech)
  results$lasso.feats$snps
  lasso.snps <- sub("_[^_]+$", "",   results$lasso.feats$snps)
  which(lasso.txpts %in% lasso.snps)
  #none overlap
  
  plot(results$lasso$X1,lgep.vs.ew.data$test.phenos$Volume)
  
  f.dat <- data.frame("Predicted"=results$lasso$X1,"True"=lgep.vs.ew.data$test.phenos$Volume)
  ggplot(f.dat,aes(x=Predicted,y=True)) + geom_point() + geom_smooth(method="lm") +
    stat_cor() + theme_pubclean()
    
    geom_vline(xintercept = mean) + geom_hline(yintercept = mean(f.dat$true.values))   
  