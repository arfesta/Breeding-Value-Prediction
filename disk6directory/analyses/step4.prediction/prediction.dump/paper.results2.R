## All files are saved as outputs and can be loaded below ####
load("/mnt/bio.012.df.RData")
load("/mnt/fam.012.df.RData")
load("/mnt/asreml.counts.RData")
load("/mnt/fam.asreml.counts.RData")
load("/mnt/bio.phenos.RData")
load("/mnt/fam.phenos.RData")
load("/mnt/fam.hm.012.df.RData")
load("/mnt/bio.hm.012.df.RData")
u.fams <- unique(fam.phenos$fam_id)
library(caret); library(caretEnsemble); library(Metrics); library(OmicKriging); library(glmnet)
# With all starting files loaded we can now generate the first plot ####
create.data.table <- function(model.list=t2) {
  
  # preallocate data types
  i = 1; MAX = length(model.list);
  x1 <- character() # Name
  x2 <- numeric()   # R2
  x3 <- numeric()   # RMSE
  x4 <- numeric()   # time [s]
  x5 <- character() # long model name
  
  # fill data and check indexes and NA
  for (i in 1:length(model.list)) {
    x1[i] <- model.list[[i]]$method
    x2[i] <- as.numeric(model.list[[i]]$results$Rsquared[which.min(model.list[[i]]$results$Rsquared)])
    x3[i] <- as.numeric(model.list[[i]]$results$RMSE[which.min(model.list[[i]]$results$RMSE)])
    x4[i] <- as.numeric(model.list[[i]]$times$everything[3])
    x5[i] <- model.list[[i]]$modelInfo$label
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
                backgroundPosition = 'center')
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
# Create even test groups ####
# Index the phenotypes so they are ind decreasing order
pheno.order <- order(fam.phenos$Volume,decreasing = T)
first <- seq(1,56,7)
last <- seq(7,56,7)

se <- c(3223,42435,23534,5235,3423,1345,21541,35432)
pheno.groups <- lapply(1:8,function(each.group){
  set.seed(se[each.group])
  sample(pheno.order[first[each.group]:last[each.group]])
})
pheno.groups <- do.call(rbind,pheno.groups)
# Each of the 7 groups contain 8 individuals spread across the population of phenotypes


### BCV ####
snp.var <- apply(fam.hm.012.df,2,var)
snp.sd <-    apply(fam.hm.012.df,2,sd)
snp.mean <-    apply(fam.hm.012.df,2,mean)

bcv.snps <- snp.sd/snp.mean
hist(bcv.snps)
# For each of the test groups, run an anova on the train.test ####
# Generate anova scores for groups 
#group.snps.anova.scores <- mclapply(1:7,function(each.group){
#  test.set <- which(bio.phenos$fam_id %in% u.fams[c(pheno.groups[,each.group])])
#  train.set <- bio.hm.012.df[-test.set,]
#  train.y <- bio.phenos$Volume[-test.set]
#  apply(train.set,2,function(each.snp) anovaScores(each.snp,train.y))
#},mc.cores=10)

#test.cts <- scale(center=F,log2(asreml.counts+1))
#group.txpt.anova.scores.v2 <- mclapply(1:7,function(each.group){
#  test.grp <- pheno.groups[,each.group]
#  test.fams <- which(bio.phenos$fam_id %in% u.fams[test.grp])
#  train.set <- test.cts[-test.fams,]
#  train.y <- bio.phenos$Volume[-test.fams]
#  apply(train.set,2,function(each.snp) anovaScores(each.snp,train.y))
#},mc.cores=20)
#group.txpt.anova.scores <- group.txpt.anova.scores.v2

#save(group.snps.anova.scores,file="/mnt/group.snps.anova.scores.RData",compress=T)
#save(group.txpt.anova.scores,file="/mnt/group.txpt.anova.scores.RData",compress=T)
load("/mnt/group.snps.anova.scores.RData")
load("/mnt/group.txpt.anova.scores.RData")
group.snps.anova.adj <- lapply(group.snps.anova.scores,function(x) p.adjust(x,method = "fdr"))
group.txpt.anova.adj <- lapply(group.txpt.anova.scores,function(x) p.adjust(x,method="fdr"))

# IMAGE-7FOLD: All txpts 7-Fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .05)]))))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(.8),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
plot(all.predictions,all.truth)
plot(group.cor); plot(group.rmse)
cor(all.predictions,all.truth)^2


gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions1 <- all.predictions
# IMAGE-7FOLD: All snps 7-fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .3)]))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
plot(all.predictions,all.truth)
plot(group.cor); plot(group.rmse)
cor(all.predictions,all.truth)^2
#gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
#gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
#all.predictions1 <- all.predictions
#IMAGE-7fOLD: All snps + All txpts with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .8)]))))
  snp.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .8)]))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.98,.02),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}

plot(all.predictions,all.truth)
cor(all.predictions,all.truth)^2

gp3.rmse <- cbind("model"=rep("transcripts*snps",7),"rmse"=group.rmse)
gp3.cor <- cbind("model"=rep("transcripts*snps",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions3 <- all.predictions
  

# Glmnet ####
pval <- .15; the.pred <- c()
for(each.group in 1:7){
test.group <- pheno.groups[,each.group]
train.group <- c(1:56)[-test.group]
snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < pval)]
#txpt.cor <- scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < pval)] + 1))
#colnames(txpt.cor) <- colnames(asreml.counts)[which(group.txpt.anova.adj[[each.group]] < pval)]

#df <- cbind(snp.cor,txpt.cor)
df.high.cor <- findCorrelation(cor(snp.cor),cutoff = .9)
snp.cor <- snp.cor[,-df.high.cor]

  set.seed(34234)
  t3 <- glmnet(y=fam.phenos$Volume[train.group], x=snp.cor[train.group,],alpha =1,standardize = F)
  prediction <- predict(t3,snp.cor[test.group,])
  the.pred <- c(the.pred,prediction[,ncol(prediction)])
  }
  
  cor(the.pred,fam.phenos$Volume[c(pheno.groups)])^2
  plot(the.pred,fam.phenos$Volume[c(pheno.groups)])
# Caret#####  
# Try varios p-thresholds for snps
{

tC <- trainControl(method="LGOCV",allowParallel = T,savePredictions="all",verboseIter = T,search="grid")
algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls","knn","leapForward","leapSeq")
pval <- .05
cv.fold.mod <- lapply(1:7,function(each.group){
test.group <- pheno.groups[,each.group]
train.group <- c(1:56)[-test.group]
snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < pval)]
df.high.cor <- findCorrelation(cor(snp.cor),cutoff = .9)
snp.cor <- snp.cor[,-df.high.cor]

t2 <- mclapply(1:length(algorithmList),function(i){ 
  set.seed(34234)
  t3 <- train(y=fam.phenos$Volume[train.group], x=snp.cor[train.group,],algorithmList[i], trControl = tC,metric = "RMSE")
},mc.cores = 20,mc.allow.recursive = T)
t2
  })

group.predictions <- lapply(1:7,function(each.group){
                                  do.call(cbind,
                                          lapply(cv.fold.mod[[each.group]], function(i)
                                            {
                                            snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < pval)]
                                            df.high.cor <- findCorrelation(cor(snp.cor))
                                            snp.cor <- snp.cor[,-df.high.cor]
                                            (predict(i,snp.cor[pheno.groups[,each.group],]))
                                            }))
  })
group.predictions <- do.call(rbind,group.predictions)
group.r2 <- apply(group.predictions,2,function(c) cor(c,fam.phenos$Volume[c(pheno.groups)]))
results.df <- as.numeric(as.character(group.predictions))
results.df<- data.frame(cbind("fam"=rep(fam.phenos$Volume[c(pheno.groups)],times=8),"prediction"=results.df))
results.df$prediction <- as.numeric(as.character(results.df$prediction))
tgc <- summarySE(results.df, measurevar="prediction",groupvars=c("fam"))
pd <- position_dodge(0.1)
ggplot(tgc, aes(x=fam, y=prediction, colour=fam, group=fam)) + 
  geom_errorbar(aes(ymin=prediction-se, ymax=prediction+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("True value") +
  ylab("Predicted Value") +
  ggtitle("Average Prediction across 12 models") +
  expand_limits(y=100) +                        # Expand y range
  scale_y_continuous() +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0))               # Position legend in bottom right 
}

# Try various p-thresholds for txpts
{
  group.txpts.anova.adj <- lapply(group.txpt.anova.scores,function(each.list) p.adjust(each.list,method = "fdr"))
  tC <- trainControl(method="LGOCV",allowParallel = T,savePredictions="all",verboseIter = T,search="grid")
  
  algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls","knn","leapForward","leapSeq")
  pval <- .05
  cv.fold.mod.txpts <-lapply(1:7,function(each.group){
    test.group <- pheno.groups[,each.group]
    train.group <- c(1:56)[-test.group]
    txpt.cor <- scale(center=F,log2(fam.asreml.counts[,which(group.txpts.anova.adj[[each.group]] < pval)] + 1))
    colnames(txpt.cor) <- colnames(asreml.counts)[which(group.txpts.anova.adj[[each.group]] < pval)]

 txpt.cor <- apply(txpt.cor,2,function(x)as.integer(x))
    t2 <- mclapply(1:length(algorithmList),function(i){ 
      set.seed(34234)
      t3 <- train(y=fam.phenos$Volume[train.group], x=txpt.cor[train.group,],algorithmList[i],trControl = tC)
    },mc.cores = 10,mc.allow.recursive = T)
    t2
  })
  
  group.predictions <- lapply(1:7,function(each.group){
    do.call(cbind,
            lapply(cv.fold.mod.txpts[[each.group]], function(i)
            {
              txpt.cor <- scale(center=F,log2(fam.asreml.counts[,which(group.txpts.anova.adj[[each.group]] < pval)] + 1))
              colnames(txpt.cor) <- colnames(asreml.counts)[which(group.txpts.anova.adj[[each.group]] < pval)]
              
              (predict(i,txpt.cor[pheno.groups[,each.group],]))}))
  })
  group.predictions <- do.call(rbind,group.predictions)
  group.r2 <- apply(group.predictions,2,function(c) cor(c,fam.phenos$Volume[c(pheno.groups)]))
  results.df <- as.numeric(as.character(group.predictions))
  results.df<- data.frame(cbind("fam"=rep(fam.phenos$Volume[c(pheno.groups)],times=8),"prediction"=results.df))
  results.df$prediction <- as.numeric(as.character(results.df$prediction))
  tgc <- summarySE(results.df, measurevar="prediction",groupvars=c("fam"))
  pd <- position_dodge(0.1)
  ggplot(tgc, aes(x=fam, y=prediction, colour=fam, group=fam)) + 
    geom_errorbar(aes(ymin=prediction-se, ymax=prediction+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("True value") +
    ylab("Predicted Value") +
    ggtitle("Average Prediction across 12 models") +
    expand_limits(y=100) +                        # Expand y range
    scale_y_continuous() +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(1,0))               # Position legend in bottom right 
}

