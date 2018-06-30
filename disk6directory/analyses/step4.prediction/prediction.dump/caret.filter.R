## SETUP ####
# Identify crossover families between batches
crossover.fams <- unique(bio.pheno.subset[,c("fam_id","Volume","batch")])
crossover.fams <- names(which(table(crossover.fams$fam_id) ==2))

# Subset the unique list of families that are in the LGEP or EW
ew.fams <- unique(bio.pheno.subset$fam_id[which(bio.pheno.subset$batch == "EW")])
lgep.fams <- unique(bio.pheno.subset$fam_id[which(bio.pheno.subset$batch == "LGEP")])

# Remove LGEP families from the EW families (i.e. remove crossover families)
ew.fams <- ew.fams[-which(ew.fams %in% crossover.fams)]

# Create a training set from the biological replicates using the lgep and ew families
train.p <- bio.pheno.subset$Volume[which(bio.pheno.subset$fam_id %in% lgep.fams)]
train.dat <- out.012.subset[which(bio.pheno.subset$fam_id %in% lgep.fams),]
colnames(train.dat) <- 1:ncol(train.dat)

# Using the biological replicate training set, apply anovaScores and gamScores to each transcript
anova_scores <- unlist(mclapply(1:ncol(train.dat),function(each.col){anovaScores(x = train.dat[,each.col],y = train.p)},mc.cores=30))
hist(anova_scores)
length(which(p.adjust(anova_scores) < .05))
#88
#gam_scores <-  unlist(mclapply(1:ncol(train.dat),function(each.col){gamScores(x = train.dat[,each.col],y = train.p)},mc.cores=30))
#hist(gam_scores)
#length(which(p.adjust(gam_scores) < .05))
#88

# * Appears that both have 88 transcripts and identify the same subset

# Setup train control and test GLMNET ####
trainControl <- trainControl(method="LOOCV",allowParallel = T, savePredictions="all",verboseIter = F,search="grid")

# For each threshold, subset the adjusted anova scores at a the given level at use for prediction of EW in caret
p.val.thresh <- seq(.01,.2,.01)

pred.results <- mclapply(1:length(p.val.thresh),function(each.thresh){
  #each.thresh=1 ## test
  train.dat.a05 <- fam.012[,which(p.adjust(anova_scores) < p.val.thresh[each.thresh])]
  train.caret.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
  train.caret.dat <- train.dat.a05[which(fam.phenos$fam_id %in% lgep.fams),]
  colnames(train.caret.dat) <- 1:ncol(train.caret.dat)
  
  mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="glmnet")
  #mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="BstLm")
  
  test.dat <- train.dat.a05; colnames(test.dat) <- colnames(train.caret.dat)
  pred <- predict(mod1,test.dat)[which(fam.phenos$fam_id %in% ew.fams)]
  
  r2 = cor(pred,fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)])^2
  rmse = rmse(actual = fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)],predicted = pred)
  
  results <- list(r2,rmse)
  return(results)
  
},mc.cores=length(p.val.thresh))

pred.r2 <- unlist(lapply(pred.results,function(x) x[1]))
pred.rmse <- unlist(lapply(pred.results,function(x) x[2]))

plot(p.val.thresh,unlist(pred.r2),xlab="p-val thresh",ylab="R2",main="Prediction of EW BV's: Adjusted p-val .01 to .2")
plot(p.val.thresh,unlist(pred.rmse),xlab="p-val thresh",ylab="RMSE",main="Prediction of EW BV's: Adjusted p-val .01 to .2")



# Setup train control and test BSTLM ####
trainControl <- trainControl(method="LOOCV",allowParallel = F, savePredictions="all",verboseIter = F,search="grid")

# For each threshold, subset the adjusted anova scores at a the given level at use for prediction of EW in caret
p.val.thresh <- seq(.01,.15,.01)

pred.results <- mclapply(1:length(p.val.thresh),function(each.thresh){
  #each.thresh=1 ## test
  train.dat.a05 <- fam.012[,which(p.adjust(anova_scores) < p.val.thresh[each.thresh])]
  train.caret.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
  train.caret.dat <- train.dat.a05[which(fam.phenos$fam_id %in% lgep.fams),]
  colnames(train.caret.dat) <- 1:ncol(train.caret.dat)
  
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="glmnet")
  mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="BstLm")
  
  test.dat <- train.dat.a05; colnames(test.dat) <- colnames(train.caret.dat)
  pred <- predict(mod1,test.dat)[which(fam.phenos$fam_id %in% ew.fams)]
  
  r2 = cor(pred,fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)])^2
  rmse = rmse(actual = fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)],predicted = pred)
  
  results <- list(r2,rmse)
  return(results)
  
},mc.cores=length(p.val.thresh))

pred.r2 <- unlist(lapply(pred.results,function(x) x[1]))
pred.rmse <- unlist(lapply(pred.results,function(x) x[2]))

plot(p.val.thresh,unlist(pred.r2),xlab="p-val thresh",ylab="R2",main="Prediction of EW BV's: Adjusted p-val .01 to .2")
plot(p.val.thresh,unlist(pred.rmse),xlab="p-val thresh",ylab="RMSE",main="Prediction of EW BV's: Adjusted p-val .01 to .2")

# Setup train control and test xgbLinear ####
trainControl <- trainControl(method="LOOCV",allowParallel = T, savePredictions="all",verboseIter = T,search="grid")

# For each threshold, subset the adjusted anova scores at a the given level at use for prediction of EW in caret
p.val.thresh <- seq(.01,.15,.01)

pred.results <- mclapply(1:length(p.val.thresh),function(each.thresh){
  each.thresh=10 ## test
  
  train.dat.a05 <- fam.012[,which(p.adjust(anova_scores) < p.val.thresh[each.thresh])]
  train.caret.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
  train.caret.dat <- train.dat.a05[which(fam.phenos$fam_id %in% lgep.fams),]
  colnames(train.caret.dat) <- 1:ncol(train.caret.dat)
  
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="glmnet")
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="BstLm")
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="glmboost")
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="svmLinear")
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="simpls")
  #mod1 <- train(y=train.caret.p,x = train.caret.dat,trControl=trainControl,method="penalized")
  test.dat <- train.dat.a05; colnames(test.dat) <- colnames(train.caret.dat)
  pred <- predict(mod1,test.dat)[which(fam.phenos$fam_id %in% ew.fams)]
  
  r2 = cor(pred,fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)])^2
  rmse = rmse(actual = fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)],predicted = pred)
  
  results <- list(r2,rmse)
  return(results)
  
},mc.cores=length(p.val.thresh))

pred.r2 <- unlist(lapply(pred.results,function(x) x[1]))
pred.rmse <- unlist(lapply(pred.results,function(x) x[2]))

plot(p.val.thresh,unlist(pred.r2),xlab="p-val thresh",ylab="R2",main="Prediction of EW BV's: Adjusted p-val .01 to .2")
plot(p.val.thresh,unlist(pred.rmse),xlab="p-val thresh",ylab="RMSE",main="Prediction of EW BV's: Adjusted p-val .01 to .2")
