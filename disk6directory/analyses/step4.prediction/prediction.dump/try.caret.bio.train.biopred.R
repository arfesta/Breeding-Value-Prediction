
LGEP.train.dat <- out.012.subset[which(bio.pheno.subset$fam_id %in% lgep.fams),which(p.adjust(anova_scores) < .05)]
colnames(LGEP.train.dat) <- colnames(out.012.subset[,which(p.adjust(anova_scores) < .05)])
EW.test.dat <- out.012.subset[which(bio.pheno.subset$fam_id %in% ew.fams),which(p.adjust(anova_scores) < .05)]
colnames(EW.test.dat) <- colnames(out.012.subset[,which(p.adjust(anova_scores) < .05)])
LGEP.train.p <- bio.pheno.subset$Volume[which(bio.pheno.subset$fam_id %in% lgep.fams)]
EW.test.p <- bio.pheno.subset$Volume[which(bio.pheno.subset$fam_id %in% ew.fams)]


# Create test/train
lgep.train <- bio.pheno.subset[which(bio.pheno.subset$fam_id %in% lgep.fams),]
test.list <- lapply(1:length(lgep.fams),function(each.lgep){
 which(lgep.train$fam_id %in% lgep.fams[each.lgep])
})

train.list <- lapply(1:length(lgep.fams),function(each.lgep){
  c(1:nrow(lgep.train))[-which(lgep.train$fam_id %in% lgep.fams[each.lgep])]
})

tC <- trainControl(method="LGOCV",allowParallel = T, savePredictions="all",verboseIter = TRUE,search="grid",index = train.list,indexOut = test.list)

library(doParallel); cl <- makeCluster(detectCores()); registerDoParallel(cl)
algorithmList <- c("glmnet","pls","bagEarth","bagEarthGCV","cforest","gaussprRadial","knn","penalized","ranger","RRFglobal","svmLinear2","svmRadialCost")

t2 <- lapply(algorithmList,function(i){ 
  set.seed(123)
 t2 <- train(y=LGEP.train.p, x=LGEP.train.dat, (i), trControl = tC)
}
)

r2 <- lapply(1:length(t2), function(i) 
{cat(sprintf("%-20s",(algorithmList[i])));
  cat(round(t2[[i]]$results$Rsquared[which.min(t2[[i]]$results$RMSE)],4),"\t");
  cat(round(t2[[i]]$results$RMSE[which.min(t2[[i]]$results$RMSE)],4),"\t")
  cat(t2[[i]]$times$everything[3],"\n")
}
)
stopCluster(cl); registerDoSEQ();


# preallocate data types
i = 1; MAX = length(t2);
x1 <- character() # Name
x2 <- numeric()   # R2
x3 <- numeric()   # RMSE
x4 <- numeric()   # time [s]
x5 <- character() # long model name

# fill data and check indexes and NA
for (i in 1:length(t2)) {
  x1[i] <- t2[[i]]$method
  x2[i] <- as.numeric(t2[[i]]$results$Rsquared[which.min(t2[[i]]$results$RMSE)])
  x3[i] <- as.numeric(t2[[i]]$results$RMSE[which.min(t2[[i]]$results$RMSE)])
  x4[i] <- as.numeric(t2[[i]]$times$everything[3])
  x5[i] <- t2[[i]]$modelInfo$label
}

# coerce to data frame
df1 <- data.frame(x1,x2,x3,x4,x5, stringsAsFactors=FALSE)

# call web browser output with sortable column names
datatable(df1,  options = list(
  columnDefs = list(list(className = 'dt-left', targets = c(0,1,2,3,4,5))),
  pageLength = MAX,
  order = list(list(2, 'desc'))),
  colnames = c('Num', 'Name', 'R^2', 'RMSE', 'time [s]', 'Model name'),
  caption = paste('Regression results from caret models',Sys.time()),
  class = 'cell-border stripe')  %>% 	       
  formatRound('x2', 3) %>%  
  formatRound('x3', 3) %>%
  formatRound('x4', 3) %>%
  formatStyle(2,
              background = styleColorBar(x2, 'steelblue'),
              backgroundSize = '100% 90%',
              backgroundRepeat = 'no-repeat',
              backgroundPosition = 'center'
  )

ew.bio.phenos <- bio.pheno.subset[which(bio.pheno.subset$fam_id %in% ew.fams),]
for(each.model in 1:length(t2)){
the.prediction <- predict(t2[[each.model]],EW.test.dat)
fam.pred <- c()
for(each.fam in 1:length(ew.fams)){
  these.rows <- bio.pheno.subset$animal_id[which(bio.pheno.subset$fam_id %in% c(ew.fams[each.fam]))]
  these.samps <- which(ew.bio.phenos$animal_id %in% these.rows)
  fam.pred <- c(fam.pred,mean(the.prediction[these.samps]))
}
print(cor(fam.pred,fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)]))
}


for(each.model in 6:length(t2)){
  the.prediction <- predict(t2[[each.model]],EW.test.dat)
  plot(the.prediction,EW.test.p)
}
predict(t2[[1]],EW.test.dat)
results.cor <- unlist(lapply(t2,function(i){cor(predict(i,EW.test.dat),EW.test.p)}))
results.rmse <- unlist(lapply(t2,function(i){RMSE(pred = predict(i,EW.test.dat),obs = EW.test.p)}))

plot(results.cor,df1$x2)
plot(results.cor,results.rmse)
this.mod <- which.min(results.rmse)
# model 3 lowest MSE
plot(predict(t2[[this.mod]],EW.test.dat),EW.test.p)

this.mod <- which.max(results.cor)
plot(predict(t2[[this.mod]],EW.test.dat),EW.test.p)
