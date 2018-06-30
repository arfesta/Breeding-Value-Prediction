install.packages(pkgs=c("caret","caretEnsemble","Metrics","DT"))

tC <- trainControl(method="LOOCV",allowParallel = F, savePredictions="all",verboseIter = F,search="grid")
alpha <- .5
  fam.preds <- mclapply(1:56,function(each.fam){
    train.dat <- fam.012[-each.fam,which( p.adjust(family_anova_scores[[each.fam]]) < alpha)]
    train.p <- fam.phenos$Volume[-each.fam]
    
    test.dat <- fam.012[,which( p.adjust(family_anova_scores[[each.fam]]) < alpha)]
    #set.seed(123)
    t2 <- train(y=train.p, x=train.dat, method = "pls", trControl = tC,scale=T)
    unlist(predict(t2,test.dat)[each.fam])
  },mc.cores=30)
  
  plot(unlist(fam.preds),fam.phenos$Volume)
  abline(h=mean(fam.phenos$Volume),v=mean(unlist(fam.preds)))
  
  
  
  library(caretEnsemble)
  # Stacking Algorithms - Run multiple algos in one call.
  alpha = .1
  algorithmList <- c('glmnet', 'svmLinear2', 'penalized','simpls')
  fam.preds <- mclapply(1:56,function(each.fam){
    train.dat <- fam.012[-each.fam,which( p.adjust(family_anova_scores[[each.fam]]) < alpha)]
    train.p <- fam.phenos$Volume[-each.fam]
    test.dat <- fam.012[,which( p.adjust(family_anova_scores[[each.fam]]) < alpha)]
    
  tc <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=2,
                               savePredictions=TRUE)
  
  set.seed(100)
  models <- caretList(y=train.p, x=train.dat, methodList = algorithmList, trControl = tc,standardize=T)  
  #results <- resamples(models)
  #summary(results)
  
  #scales <- list(x=list(relation="free"), y=list(relation="free"))
  #bwplot(results, scales=scales)
  
  set.seed(101)
  stackControl <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=2,
                               savePredictions=TRUE)
  
  # Ensemble the predictions of `models` to form a new combined prediction based on glm
  stack.glm <- caretStack(models, method="glm", metric="RMSE", trControl=stackControl)
  predict(stack.glm, newdata=test.dat)[each.fam]},mc.cores=30)
  
  plot(unlist(fam.preds),fam.phenos$Volume)
  