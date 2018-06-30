# Results for Paper
## All files are saved as outputs and can be loaded below ####
#load("/mnt/bio.012.df.RData")
#load("/mnt/fam.012.df.RData")
#load("/mnt/asreml.counts.RData")
load("/mnt/fam.asreml.counts.RData")
#load("/mnt/bio.phenos.RData")
load("/mnt/fam.phenos.RData")
load("/mnt/fam.hm.012.df.RData")
#load("/mnt/bio.hm.012.df.RData")
u.fams <- unique(fam.phenos$fam_id)
load("/mnt/group.snps.anova.scores.nolowvar.RData")
load("/mnt/group.txpt.anova.scores.nolowvar.RData")
group.snps.anova.adj <- lapply(group.snps.anova.scores,function(x) p.adjust(x,method = "fdr"))
group.txpt.anova.adj <- lapply(group.txpt.anova.scores,function(x) p.adjust(x,method="fdr"))
library(caret); library(caretEnsemble); library(Metrics); library(OmicKriging); library(glmnet)

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



# Caret#####  
# Try varios p-thresholds for snps
{
  
  tC <- trainControl(method="LGOCV",allowParallel = T,savePredictions="all",verboseIter = T,search="grid")
  algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls","knn","leapForward","leapSeq")
  pval.list <- seq(.01,.2,.01); output.list <- vector("list")
  for(each.pval in 1:length(pval.list)){
  pval <- pval.list[each.pval]
  cv.fold.mod <- mclapply(1:7,function(each.group){
    test.group <- pheno.groups[,each.group]
    train.group <- c(1:56)[-test.group]
    snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < pval)]
    df.high.cor <- findCorrelation(cor(snp.cor),cutoff = .9)
    snp.cor <- snp.cor[,-df.high.cor]
    t2 <- mclapply(1:length(algorithmList),function(i){ 
      set.seed(34234)
      t3 <- train(y=fam.phenos$Volume[train.group], x=snp.cor[train.group,],algorithmList[i], trControl = tC,metric = "RMSE")
    },mc.cores = 15,mc.allow.recursive = T)
    t2
  },mc.cores=2,mc.allow.recursive = T)
  output.list[[each.pval]] <- cv.fold.mod
  print(pval.list[each.pval])
  }
  
  save(output.list,file = "/mnt/snps.8mods.20pvals.RData",compress=T)
  
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
  tgc <- summarySE(all.dat.7fold, measurevar="cor",groupvars=c("model","group"))
  pd <- position_dodge(0.1)
  ggplot(tgc, aes(x=model, y=cor, colour=model)) + 
    geom_errorbar(aes(ymin=cor-se, ymax=cor+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("True value") +
    ylab("Predicted Value") +
    ggtitle("Average Prediction across 12 models") +
    expand_limits() +                        # Expand y range
    scale_y_continuous() +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(1,0))               # Position legend in bottom right 
}


{
colnames(fam.asreml.counts) <- colnames(asreml.counts)
  test.cts <- scale(center=F,log2(fam.asreml.counts+1))
  tC <- trainControl(method="LGOCV",allowParallel = T,savePredictions="all",verboseIter = T,search="grid")
  algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls","knn","leapForward","leapSeq")
  pval.list <- seq(.01,.2,.01); output.list.txpts <- vector("list")
  for(each.pval in 1:length(pval.list)){
    pval <- pval.list[each.pval]
    cv.fold.mod <- mclapply(1:7,function(each.group){
      test.group <- pheno.groups[,each.group]
      train.group <- c(1:56)[-test.group]
      txpt.cor <- test.cts[,which(group.txpt.anova.adj[[each.group]] < pval)]
      df.high.cor <- findCorrelation(cor(txpt.cor),cutoff = .9)
      txpt.cor <- txpt.cor[,-df.high.cor]
      t2 <- mclapply(1:length(algorithmList),function(i){ 
        set.seed(34234)
        t3 <- train(y=fam.phenos$Volume[train.group], x=txpt.cor[train.group,],algorithmList[i], trControl = tC,metric = "RMSE")
      },mc.cores = 15,mc.allow.recursive = T)
      t2
    },mc.cores=2,mc.allow.recursive = T)
      output.list.txpts[[each.pval]] <- cv.fold.mod
    print(pval.list[each.pval])
  }
  
  save(output.list,file = "/mnt/txpts.8mods.20pvals.RData",compress=T)
  
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

# ALL TRANSCRIPTS OR SNPS ####
library(OmicKriging); library(Metrics)
txpt.cor <- cor(t(scale(log2(fam.asreml.counts+1),center=F)))
rownames(txpt.cor) <- rownames(fam.phenos); colnames(txpt.cor) <- rownames(fam.phenos)

snp.cor <- cor(t(fam.hm.012.df))
#snp.cor <- cor(t(fam.012.df))
rownames(snp.cor) <- rownames(fam.phenos); colnames(snp.cor) <- rownames(fam.phenos)

# Identify crossover families between batches
crossover.fams <- unique(bio.phenos[,c("fam_id","Volume","batch")])
crossover.fams <- names(which(table(crossover.fams$fam_id) ==2))
# Subset the unique list of families that are in the LGEP or EW
ew.fams <- unique(bio.phenos$fam_id[which(bio.phenos$batch == "EW")])
lgep.fams <- unique(bio.phenos$fam_id[which(bio.phenos$batch == "LGEP")])
# Remove LGEP families from the EW families (i.e. remove crossover families)
ew.fams <- ew.fams[-which(ew.fams %in% crossover.fams)]

# Set up random dataframe 
set.seed(424)
random.txpt.df <- fam.asreml.counts[sample(1:56,56),]
rownames(random.txpt.df) <- fam.phenos$fam_id
txpt.random.cor <- cor(t(scale(log2(random.txpt.df + 1),center=F)))
colnames(txpt.random.cor) <- fam.phenos$fam_id; rownames(txpt.random.cor) <- fam.phenos$fam_id

set.seed(424)
random.snp.df <- fam.012.df[sample(1:56,56),]
rownames(random.snp.df) <- fam.phenos$fam_id
snp.random.cor <- cor(t(random.snp.df))
colnames(snp.random.cor) <- fam.phenos$fam_id; rownames(snp.random.cor) <- fam.phenos$fam_id


# IMAGE-7FOLD: All txpts 7-Fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
plot(all.predictions,all.truth)
plot(group.cor); 
#plot(group.rmse)
cor(all.predictions,all.truth)^2
gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions1 <- all.predictions
# IMAGE-7FOLD: All snps 7-fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]-mean(fam.phenos$Volume[-c(pheno.groups[,each.group])])); all.truth <- c(all.truth,pred[,3]-mean(fam.phenos$Volume[-c(pheno.groups[,each.group])]))
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,3],pred[,2]))
  plot(pred[,2],pred[,3])
}
plot(all.predictions,all.truth)
plot(group.cor); 
#plot(group.rmse)
cor(all.predictions,all.truth)^2
gp2.rmse <- cbind("model"=rep("snps",7),"rmse"=group.rmse)
gp2.cor <- cbind("model"=rep("snps",7),"cor"=group.cor,"rmse"=group.rmse)   
all.predictions2 <- all.predictions
# IMAGE-7fOLD: All snps + All txpts with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.98,.02),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}

plot(all.predictions,all.truth)
gp3.rmse <- cbind("model"=rep("transcripts*snps",7),"rmse"=group.rmse)
gp3.cor <- cbind("model"=rep("transcripts*snps",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions3 <- all.predictions

## 7-FOLD Image: ####
rm.mods <- c(3,4,6,7,8)
pval.thresh <- 3
#cor.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.cor[,-rm.mods]
#spear.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.spear[,-rm.mods]
#rmse.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.rmse[,-rm.mods]
cor.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.cor
spear.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.spear
rmse.dat <- predictions.extracted.1.20[[pval.thresh]]$stats_by_group$group.rmse
all.dat.7fold <- data.frame(
  "cor" = c(cor.dat),
  "spear" = c(spear.dat),
  "rmse"= c(rmse.dat),
  "model" =as.factor(as.character(rep(x = 1:ncol(cor.dat),each=7))),
  "group" = as.factor(as.character(rep(x = 1:7,times=ncol(cor.dat))))
,stringsAsFactors = F)

qplot(model, cor, data = all.dat.7fold, 
      geom= "violin", fill = model)
qplot(model, rmse, data = all.dat.7fold, 
      geom= "violin", fill = model)

ggboxplot(all.dat.7fold, x = "group", y = "rmse",
          color = "group", palette = "jco",short.panel.labs = F,show.legend=F) + geom_point(aes(shape=model),size=2) 

#ggplot(data=all.dat.7fold,aes(x=group, y=cor,colour=model)) + geom_point(size=3) + 
#  geom_path(data = all.dat.7fold, aes(x=as.numeric(as.character(group)),y=cor)) 

#ggboxplot(all.dat.7fold, x = "group", y = "cor",
#          color = "model", palette = "jco",
#          add = "jitter",short.panel.labs = F,shape="model",size=3)
#p + ggplot(data=all.dat.7fold,aes(x=group, y=rmse,colour=model)) + geom_point(size=3) + 
#  geom_path(data = all.dat.7fold, aes(x=as.numeric(as.character(group)),y=rmse))

#final.dat <- data.frame(cbind("true.val"=rep(all.truth,times=3), 
#                              "pred.val"=c(all.predictions1,all.predictions2,all.predictions3),
#                              "model"=c(rep("transcripts",56),rep("snps",56),rep("transcripts*snps",56))))
#final.dat$true.val <- as.numeric(as.character(final.dat$true.val))
#final.dat$pred.val <- as.numeric(as.character(final.dat$pred.val))
#(plot5 <- ggplot(aes(x = true.val, y = pred.val, group = model, colour = model),
#                 data = final.dat) +
#    geom_line() + geom_point())
#t.test(as.numeric(gp1.cor[,2]),as.numeric(gp2.cor[,2]),paired = T)
#plot(all.dat.7fold$cor,all.dat.7fold$rmse)

# IMAGE-7FOLD: All txpts 7-Fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){

 txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj2[[each.group]] < .15)]))))

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


gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions1 <- all.predictions
# IMAGE-7FOLD: All snps 7-fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .05)]))
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
  txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .15)]))))
  snp.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .15)]))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.8,.2),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}

plot(all.predictions,all.truth)
cor(all.predictions,all.truth)^2

gp3.rmse <- cbind("model"=rep("transcripts*snps",7),"rmse"=group.rmse)
gp3.cor <- cbind("model"=rep("transcripts*snps",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions3 <- all.predictions



tC <- trainControl(method="LGOCV",allowParallel = T,p=.75,savePredictions="all",verboseIter = F,search="grid")

algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls")

the.predictions <- lapply(1:7,function(each.group){
  #txpt.cor <- scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .05)]))
  snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .1)]
  test.group <- c(pheno.groups[,each.group])
  train.group <- c(1:56)[-pheno.groups[,each.group]]
  t3 <- train(y=fam.phenos$Volume[train.group], x=snp.cor[train.group,],algorithmList[each.group],metric = "RMSE", trControl = tC)
  cor(predict(object = t3,newdata = snp.cor[test.group,]),fam.phenos$Volume[test.group])
  
  
  
  
  t2 <- mclapply(1:length(algorithmList),function(i){ 
    t3 <- train(y=LGEP.train.p, x=LGEP.train.dat,algorithmList[i], trControl = tC)
  },mc.cores=20)
  
  
  